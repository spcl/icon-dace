/**
 * @file ppm_checksum_c.c
 * @brief support checksumming of data in C
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>

#ifdef USE_CRYPTO
#include <openssl/evp.h>
#include <openssl/opensslv.h>
#endif

#include "core/core.h"
#include "crypto/ppm_checksum.h"

#include "core/ppm_visibility.h"
#ifndef USE_CRYPTO
#include "crypto/md5.h"
#endif

void PPM_DSO_API_EXPORT
PPM_describe_digest(enum digest_type digest, struct PPM_digest_description *hd)
{
  switch (digest)
  {
#ifdef USE_CRYPTO
  case digest_type_md5:
    hd->helper = (void *)EVP_md5();
    hd->size = (int32_t)EVP_MD_size((const EVP_MD *)(hd->helper));
    hd->family = digest_crypto;
    break;
  case digest_type_sha1:
    hd->helper = (void *)EVP_sha1();
    hd->size = (int32_t)EVP_MD_size((const EVP_MD *)(hd->helper));
    hd->family = digest_crypto;
    break;
#else
  case digest_type_md5:
    hd->helper = (void *)(intptr_t)digest_type_md5;
    hd->size = PPM_MD5_DIGEST_LENGTH;
    hd->family = digest_internal;
    break;
  case digest_type_sha1:
    hd->helper = NULL;
    hd->size = 0;
    hd->family = digest_none;
    break;
#endif
  default:
    PPM_abort(PPM_default_comm, "invalid digest type specified",
              __FILE__, __LINE__);
  }
}

void PPM_DSO_API_EXPORT
PPM_checksum(const void *buf, size_t buf_size, unsigned char *checksum,
             struct PPM_digest_description *digest)
{
  assert(((buf_size > 0 && buf ) || !buf_size) && checksum);
  switch (digest->family) {
  case digest_none:
    PPM_abort(PPM_default_comm, "invalid digest type specified",
              __FILE__, __LINE__);
    break;
  case digest_crypto:
#ifdef USE_CRYPTO
    {
      unsigned int md_len;
      const EVP_MD *md = digest->helper;
#if OPENSSL_VERSION_NUMBER < 0x1010000fL
      EVP_MD_CTX mdctx;
      EVP_MD_CTX_init(&mdctx);
#define pmdctx &mdctx
#define EVP_MD_CTX_free(ctx)
#else
      EVP_MD_CTX *pmdctx = EVP_MD_CTX_new();
      if (pmdctx == NULL) {
	  PPM_abort(PPM_default_comm, "failed to initialize digest context",
		    __FILE__, __LINE__);
      }
#endif
      if (!EVP_DigestInit_ex(pmdctx, md, NULL)) {
	  EVP_MD_CTX_free(pmdctx);
	  PPM_abort(PPM_default_comm, "failed to setup digest context",
		    __FILE__, __LINE__);
      }
      if (!EVP_DigestUpdate(pmdctx, buf, buf_size)) {
	  EVP_MD_CTX_free(pmdctx);
	  PPM_abort(PPM_default_comm, "failed to hash data into digest context",
		    __FILE__, __LINE__);
      }
      if (!EVP_DigestFinal_ex(pmdctx, checksum, &md_len)) {
	  EVP_MD_CTX_free(pmdctx);
	  PPM_abort(PPM_default_comm, "failed to retrieve digest value",
		    __FILE__, __LINE__);
      }
      EVP_MD_CTX_free(pmdctx);
    }
#else
    PPM_abort(PPM_default_comm, "invalid digest type specified",
              __FILE__, __LINE__);
#endif
    break;
  case digest_internal:
#ifndef USE_CRYPTO
    switch ((int)(intptr_t)digest->helper) {
    case digest_type_md5:
      {
        PPM_MD5_CTX mdctx;
        PPM_MD5Init(&mdctx);
        PPM_MD5Update(&mdctx, buf, buf_size);
        PPM_MD5Final(checksum, &mdctx);
      }
      return;
    }
#endif
  default:
    PPM_abort(PPM_default_comm, "invalid digest type specified",
              __FILE__, __LINE__);
  }
}

void PPM_DSO_API_EXPORT
PPM_md2hex(char *hex, const unsigned char *checksum, size_t checksum_size)
{
  size_t i;
  for (i = 0; i < checksum_size; ++i)
    sprintf(hex + 2 * i, "%02x", checksum[i]);
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
