#include "cdi_across.h"

#ifdef HAVE_ACROSS

#include <assert.h>
#include <netdb.h>
#include <string.h>
#include <sys/socket.h>
#include <time.h>
#include <unistd.h>

#include "binary.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "error.h"

#define ACROSS_SECTION2_VERSION 1
#define EXPID_MIN_LENGTH 7

static uint16_t
swap_endianness16(uint16_t v)
{
  // clang-format off
  return ((v >> 8) & 0xff)
       | ((v << 8) & 0xff00);
  // clang-format on
}

static uint32_t
swap_endianness32(uint32_t v)
{
  // clang-format off
  return ((v >> 24) & 0xff)
       | ((v <<  8) & 0xff0000)
       | ((v >>  8) & 0xff00)
       | ((v << 24) & 0xff000000);
  // clang-format on
}

static uint64_t
swap_endianness64(uint64_t v)
{
  // clang-format off
  return ((v >> 56) & 0xff)
       | ((v >> 40) & 0xff00)
       | ((v >> 24) & 0xff0000)
       | ((v >>  8) & 0xff000000)
       | ((v <<  8) & 0xff00000000)
       | ((v << 24) & 0xff0000000000)
       | ((v << 40) & 0xff000000000000)
       | ((v << 56) & 0xff00000000000000);
  // clang-format on
}

static int
set8(unsigned char *buf, uint8_t v)
{
  *buf = v;
  return sizeof(v);
}

static int
set16(unsigned char *buf, uint16_t v)
{
  if (HOST_ENDIANNESS == CDI_LITTLEENDIAN) v = swap_endianness16(v);
  memcpy(buf, &v, sizeof(v));
  return sizeof(v);
}

static int
set32(unsigned char *buf, uint32_t v)
{
  if (HOST_ENDIANNESS == CDI_LITTLEENDIAN) v = swap_endianness32(v);
  memcpy(buf, &v, sizeof(v));
  return sizeof(v);
}

static int
set64(unsigned char *buf, uint64_t v)
{
  if (HOST_ENDIANNESS == CDI_LITTLEENDIAN) v = swap_endianness64(v);
  memcpy(buf, &v, sizeof(v));
  return sizeof(v);
}

static void
fletcher8(const unsigned char *buf, size_t len, uint8_t *sum1out, uint8_t *sum2out)
{
  uint16_t sum1 = 0;
  uint16_t sum2 = 0;
  for (size_t i = 0; i < len; ++i)
    {
      sum1 += *buf++;
      if (sum1 >= 255) sum1 -= 255;
      sum2 += sum1;
      if (sum2 >= 255) sum2 -= 255;
    }
  *sum1out = sum1;
  *sum2out = sum2;
}

static uint8_t
fletcher8_check1(uint8_t sum1, uint8_t sum2)
{
  return 255 - ((sum1 + sum2) % 255);
}

static uint8_t
fletcher8_check2(uint8_t sum1, uint8_t sum2)
{
  return 255 - ((sum1 + fletcher8_check1(sum1, sum2)) % 255);
}

static int
across_write_buf(int sock, const void *buf, size_t nbytes)
{
  const char *bufpos = (const char *) buf;
  size_t nbytes_left = nbytes;
  while (nbytes_left > 0)
    {
      ssize_t nbytes_written = write(sock, bufpos, nbytes_left);
      if (nbytes_written < 0)
        {
          perror(__func__);
          return 1;
        }
      nbytes_left -= nbytes_written;
      bufpos += nbytes_written;
    }
  return 0;
}

int
across_write_grib_message(stream_t *streamptr, const void *gribbuffer, size_t nbytes)
{
  unsigned char section0[] = { 'G', 'R', 'I', 'B', 255, 255, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0 };
  const size_t section0_len = sizeof(section0);

  uint32_t section1_len;
  const unsigned char *section1 = (unsigned char *) gribbuffer + section0_len;
  if (nbytes < section0_len + sizeof(section1_len) + 1)
    {
      Error("GRIB2 stream is too short");
      return 1;
    }
  if (section1[sizeof(section1_len)] != 1)
    {
      Error("Section 1 not found at expected position in GRIB2 stream");
      return 1;
    }
  memcpy(&section1_len, section1, sizeof(section1_len));
  if (HOST_ENDIANNESS == CDI_LITTLEENDIAN) section1_len = swap_endianness32(section1_len);

  if (nbytes < section0_len + section1_len + 4  // next section size
                   + 1                          // section number ("3")
                   + 4                          // footer "7777"
  )
    {
      Error("GRIB2 stream is too short");
      return 1;
    }
  size_t other_sections_len = nbytes - section0_len - section1_len;
  const unsigned char *other_sections = section1 + section1_len;
  if (other_sections[4] != 3)
    {
      if (other_sections[4] == 2)
        Error("Section 2 must not already be present in GRIB2 stream");
      else
        Error("Section 3 not found at expected position in GRIB2 stream");
      return 1;
    }

  across_info_t *info = (across_info_t *) streamptr->protocolData;

  int expid_len = strlen(info->expid);
  if (expid_len < EXPID_MIN_LENGTH)
    {
      Error("expid is too short");
      return 1;
    }
  if (expid_len > 255)
    {
      Error("expid is too long");
      return 1;
    }
  uint32_t section2_len = 4                // section size
                          + 1              // section number ("2")
                          + 2              // local section number ("0")
                          + 1              // across local section version
                          + 1              // expid length
                          + expid_len + 4  // expver
                          + 4              // timestamp
                          + 2              // check bytes
      ;
  unsigned char *section2 = (unsigned char *) Malloc(section2_len);
  {
    unsigned char pos = 0;
    pos += set32(section2 + pos, section2_len);
    pos += set8(section2 + pos, 2);   // section number
    pos += set16(section2 + pos, 0);  // local section number
    pos += set8(section2 + pos, ACROSS_SECTION2_VERSION);
    pos += set8(section2 + pos, expid_len);
    memcpy(section2 + pos, info->expid, expid_len);
    pos += expid_len;
    pos += set32(section2 + pos, info->expver);
    pos += set32(section2 + pos, time(NULL));
    uint8_t fletchersum1;
    uint8_t fletchersum2;
    fletcher8(section2 + 5, section2_len - 7, &fletchersum1, &fletchersum2);
    pos += set8(section2 + pos, fletcher8_check1(fletchersum1, fletchersum2));
    pos += set8(section2 + pos, fletcher8_check2(fletchersum1, fletchersum2));
    assert(pos == section2_len);
  }

  set64(section0 + section0_len - 8, nbytes + section2_len);  // total message length

  int res = 0;
  if (across_write_buf(streamptr->fileID, section0, section0_len) || across_write_buf(streamptr->fileID, section1, section1_len)
      || across_write_buf(streamptr->fileID, section2, section2_len)
      || across_write_buf(streamptr->fileID, other_sections, other_sections_len))
    res = 1;

  Free(section2);
  return res;
}

int
across_connect_int(char *path, char filemode, stream_t *streamptr)
{
  // parse ACROSS address format: across://<host>[:<port>]/<expid>[@<expver>]

  const char *host = path;
  const char *port = ACROSS_DEFAULT_PORT;
  const char *expid = NULL;
  int expver = 1;

  char *c = path;
  while (*c != '\0')
    {
      switch (*c)
        {
        case ':':
          if (expid == NULL)
            {
              port = c + 1;
              *c = '\0';
            }
          ++c;
          break;

        case '/':
          if (expid != NULL)
            {
              Warning("Experiment id must not contain '/'");
              return CDI_EINVAL;
            }
          expid = c + 1;
          *c = '\0';
          ++c;
          break;

        case '@':
          if (expid == NULL)
            {
              ++c;
              break;
            }
          *c = '\0';
          expver = strtoul(c + 1, &c, 10);
          if (*c != '\0')
            {
              Warning("Experiment version must be an integer");
              return CDI_EINVAL;
            }
          break;

        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
        case 'a':
        case 'b':
        case 'c':
        case 'd':
        case 'e':
        case 'f':
        case 'g':
        case 'h':
        case 'i':
        case 'j':
        case 'k':
        case 'l':
        case 'm':
        case 'n':
        case 'o':
        case 'p':
        case 'q':
        case 'r':
        case 's':
        case 't':
        case 'u':
        case 'v':
        case 'w':
        case 'x':
        case 'y':
        case 'z':
        case 'A':
        case 'B':
        case 'C':
        case 'D':
        case 'E':
        case 'F':
        case 'G':
        case 'H':
        case 'I':
        case 'J':
        case 'K':
        case 'L':
        case 'M':
        case 'N':
        case 'O':
        case 'P':
        case 'Q':
        case 'R':
        case 'S':
        case 'T':
        case 'U':
        case 'V':
        case 'W':
        case 'X':
        case 'Y':
        case 'Z':
        case '-':
        case '.':
        case '_':
        case '~': ++c; break;

        default:
          if (expid != NULL)
            {
              Warning("Experiment id must not contain '%c'", c);
              return CDI_EINVAL;
            }
          ++c;
          break;
        }
    }

  if (expid == NULL || expid[0] == '\0')
    {
      Warning("No experiment id given");
      return CDI_EINVAL;
    }

  if (strlen(expid) < EXPID_MIN_LENGTH)
    {
      Warning("Experiment id must be longer than %d chars", EXPID_MIN_LENGTH - 1);
      return CDI_EINVAL;
    }
  if (strlen(expid) > 255)
    {
      Warning("Experiment id must be shorter than 256 chars");
      return CDI_EINVAL;
    }

  struct addrinfo hints;
  memset(&hints, 0, sizeof(hints));
  hints.ai_family = AF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_flags = 0;
  hints.ai_protocol = 0;

  struct addrinfo *addr;
  int res = getaddrinfo(host, port, &hints, &addr);
  if (res != 0)
    {
      Warning("Could not find %s: %s", host, gai_strerror(res));
      return CDI_ESYSTEM;
    }

  int sock;
  struct addrinfo *rp;
  for (rp = addr; rp != NULL; rp = rp->ai_next)
    {
      sock = socket(rp->ai_family, rp->ai_socktype, rp->ai_protocol);
      if (sock < 0) continue;
      if (connect(sock, rp->ai_addr, rp->ai_addrlen) != -1) break;
      close(sock);
    }
  freeaddrinfo(addr);

  if (rp == NULL)
    {
      Warning("Could not connect to %s:%s", host, port);
      return CDI_ESYSTEM;
    }

  across_info_t *info = (across_info_t *) Malloc(sizeof(across_info_t));
  info->expid = strdup(expid);
  info->expver = expver;
  streamptr->protocolData = info;

  return sock;
}

int
across_connect(const char *path, char filemode, stream_t *streamptr)
{
  if (filemode != 'w')
    {
      Warning("Reading from ACROSS not implemented yet");  // TODO read from ACROSS
      return CDI_EINVAL;
    }

  char *p = strdup(path);  // needed to adjust part strings
  int sock = across_connect_int(p, filemode, streamptr);
  free(p);
  return sock;
}

void
across_disconnect(int sock)
{
  close(sock);
}

#endif

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
