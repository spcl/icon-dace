#ifndef __DACE_SERDE__
#define __DACE_SERDE__

#include <cassert>
#include <istream>
#include <iostream>
#include <sstream>
#include <optional>

#include "add_aerosol_optics.h"

namespace serde {
    struct array_meta {
      int rank = 0;
      std::vector<int> size, lbound;

      int volume() const {  return std::reduce(size.begin(), size.end(), 1, std::multiplies<int>()) ; }
    };
    std::map<void*, array_meta>& ARRAY_META_DICT() {
        static auto* M = new std::map<void*, array_meta>();
        return *M;
    }

    std::string scroll_space(std::istream& s) {
        std::string out;
        while (!s.eof() && (!s.peek() || isspace(s.peek()))) {
            out += s.get();
            assert(s.good());
        }
        return out;
    }

    std::string read_line(std::istream& s, const std::optional<std::string>& should_contain = {}) {
        if (s.eof()) return "<eof>";
        scroll_space(s);
        char bin[101];
        s.getline(bin, 100);
        assert(s.good());
        if (should_contain) {
            bool ok = (std::string(bin).find(*should_contain) != std::string::npos);
            if (!ok) {
                std::cerr << "Expected: '" << *should_contain << "'; got: '" << bin << "'" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        return {bin};
    }

    template<typename T>
    void read_scalar(T& x, std::istream& s) {
        if (s.eof()) return;
        scroll_space(s);
        s >> x;
    }

    void read_scalar(float& x, std::istream& s) {
        if (s.eof()) return;
        scroll_space(s);
        long double y;
        s >> y;
        x = y;
    }

    void read_scalar(double& x, std::istream& s) {
        if (s.eof()) return;
        scroll_space(s);
        long double y;
        s >> y;
        x = y;
    }

    void read_scalar(bool& x, std::istream& s) {
        char c;
        read_scalar(c, s);
        assert (c == '1' or c == '0');
        x = (c == '1');
    }

    array_meta read_array_meta(std::istream& s) {
        array_meta m;
        read_line(s, {"# rank"});  // Should contain '# rank'
        read_scalar(m.rank, s);
        m.size.resize(m.rank);
        m.lbound.resize(m.rank);
        read_line(s, {"# size"});  // Should contain '# size'
        for (int i=0; i<m.rank; ++i) {
            read_scalar(m.size[i], s);
        }
        read_line(s, {"# lbound"});  // Should contain '# lbound'
        for (int i=0; i<m.rank; ++i) {
            read_scalar(m.lbound[i], s);
        }
        return m;
    }

    
void deserialize(float* x, std::istream& s) {
    read_scalar(*x, s);
}
void deserialize(double* x, std::istream& s) {
    read_scalar(*x, s);
}
void deserialize(long double* x, std::istream& s) {
    read_scalar(*x, s);
}
void deserialize(int* x, std::istream& s) {
    read_scalar(*x, s);
}
void deserialize(long* x, std::istream& s) {
    read_scalar(*x, s);
}
void deserialize(long long* x, std::istream& s) {
    read_scalar(*x, s);
}
void deserialize(bool* x, std::istream& s) {
    read_scalar(*x, s);
}


void deserialize(aerosol_type* x, std::istream& s) {
    bool yep;
    array_meta m;
    read_line(s, {"# od_sw"});  // Should contain '# od_sw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_od_sw_d_0_s_19 = m.size[0];
x->__f2dace_SA_od_sw_d_1_s_20 = m.size[1];
x->__f2dace_SA_od_sw_d_2_s_21 = m.size[2];
x->__f2dace_SOA_od_sw_d_0_s_19 = m.lbound[0];
x->__f2dace_SOA_od_sw_d_1_s_20 = m.lbound[1];
x->__f2dace_SOA_od_sw_d_2_s_21 = m.lbound[2];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->od_sw = new std::remove_pointer<decltype(x ->od_sw)>::type[m.volume()];
ARRAY_META_DICT()[x->od_sw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->od_sw[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# ssa_sw"});  // Should contain '# ssa_sw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_ssa_sw_d_0_s_22 = m.size[0];
x->__f2dace_SA_ssa_sw_d_1_s_23 = m.size[1];
x->__f2dace_SA_ssa_sw_d_2_s_24 = m.size[2];
x->__f2dace_SOA_ssa_sw_d_0_s_22 = m.lbound[0];
x->__f2dace_SOA_ssa_sw_d_1_s_23 = m.lbound[1];
x->__f2dace_SOA_ssa_sw_d_2_s_24 = m.lbound[2];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->ssa_sw = new std::remove_pointer<decltype(x ->ssa_sw)>::type[m.volume()];
ARRAY_META_DICT()[x->ssa_sw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->ssa_sw[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# g_sw"});  // Should contain '# g_sw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_g_sw_d_0_s_25 = m.size[0];
x->__f2dace_SA_g_sw_d_1_s_26 = m.size[1];
x->__f2dace_SA_g_sw_d_2_s_27 = m.size[2];
x->__f2dace_SOA_g_sw_d_0_s_25 = m.lbound[0];
x->__f2dace_SOA_g_sw_d_1_s_26 = m.lbound[1];
x->__f2dace_SOA_g_sw_d_2_s_27 = m.lbound[2];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->g_sw = new std::remove_pointer<decltype(x ->g_sw)>::type[m.volume()];
ARRAY_META_DICT()[x->g_sw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->g_sw[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# od_lw"});  // Should contain '# od_lw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_od_lw_d_0_s_28 = m.size[0];
x->__f2dace_SA_od_lw_d_1_s_29 = m.size[1];
x->__f2dace_SA_od_lw_d_2_s_30 = m.size[2];
x->__f2dace_SOA_od_lw_d_0_s_28 = m.lbound[0];
x->__f2dace_SOA_od_lw_d_1_s_29 = m.lbound[1];
x->__f2dace_SOA_od_lw_d_2_s_30 = m.lbound[2];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->od_lw = new std::remove_pointer<decltype(x ->od_lw)>::type[m.volume()];
ARRAY_META_DICT()[x->od_lw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->od_lw[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# ssa_lw"});  // Should contain '# ssa_lw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_ssa_lw_d_0_s_31 = m.size[0];
x->__f2dace_SA_ssa_lw_d_1_s_32 = m.size[1];
x->__f2dace_SA_ssa_lw_d_2_s_33 = m.size[2];
x->__f2dace_SOA_ssa_lw_d_0_s_31 = m.lbound[0];
x->__f2dace_SOA_ssa_lw_d_1_s_32 = m.lbound[1];
x->__f2dace_SOA_ssa_lw_d_2_s_33 = m.lbound[2];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->ssa_lw = new std::remove_pointer<decltype(x ->ssa_lw)>::type[m.volume()];
ARRAY_META_DICT()[x->ssa_lw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->ssa_lw[i]), s);
}

}  // CONCLUDING IF
}


void deserialize(gas_type* x, std::istream& s) {
    bool yep;
    array_meta m;
    
}


void deserialize(config_type* x, std::istream& s) {
    bool yep;
    array_meta m;
    read_line(s, {"# i_band_from_reordered_g_lw"});  // Should contain '# i_band_from_reordered_g_lw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_i_band_from_reordered_g_lw_d_0_s_2 = m.size[0];
x->__f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_2 = m.lbound[0];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->i_band_from_reordered_g_lw = new std::remove_pointer<decltype(x ->i_band_from_reordered_g_lw)>::type[m.volume()];
ARRAY_META_DICT()[x->i_band_from_reordered_g_lw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->i_band_from_reordered_g_lw[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# i_band_from_reordered_g_sw"});  // Should contain '# i_band_from_reordered_g_sw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_i_band_from_reordered_g_sw_d_0_s_3 = m.size[0];
x->__f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_3 = m.lbound[0];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->i_band_from_reordered_g_sw = new std::remove_pointer<decltype(x ->i_band_from_reordered_g_sw)>::type[m.volume()];
ARRAY_META_DICT()[x->i_band_from_reordered_g_sw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->i_band_from_reordered_g_sw[i]), s);
}

}  // CONCLUDING IF
}


void deserialize(thermodynamics_type* x, std::istream& s) {
    bool yep;
    array_meta m;
    
}

    
template<typename T>
void add_line(const T& x, std::ostream& s, bool trailing_newline=true) {
    s << x;
    if (trailing_newline) s << std::endl;
}
void add_line(long long x, std::ostream& s, bool trailing_newline=true) {
    s << x;
    if (trailing_newline) s << std::endl;
}
void add_line(long double x, std::ostream& s, bool trailing_newline=true) {
    s << x;
    if (trailing_newline) s << std::endl;
}
void add_line(bool x, std::ostream& s, bool trailing_newline=true) {
    add_line(int(x), s, trailing_newline);
}
template<typename T>
std::string serialize(const T* x) {
    std::stringstream s;
    add_line(*x, s, false);
    return s.str();
}
std::string serialize(int x) {
    return std::to_string(x);
}
std::string serialize(long x) {
    return std::to_string(x);
}
std::string serialize(long long x) {
    return std::to_string(x);
}
std::string serialize(float x) {
    return std::to_string(x);
}
std::string serialize(double x) {
    return std::to_string(x);
}
std::string serialize(long double x) {
    return std::to_string(x);
}
std::string serialize(bool x) {
    return serialize(int(x));
}


std::string serialize(const aerosol_type* x) {
    std::stringstream s;
    add_line("# od_sw", s);

add_line("# alloc", s);
add_line(serialize(x->od_sw != nullptr), s);
if (x->od_sw) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->od_sw];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->od_sw[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# ssa_sw", s);

add_line("# alloc", s);
add_line(serialize(x->ssa_sw != nullptr), s);
if (x->ssa_sw) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->ssa_sw];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->ssa_sw[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# g_sw", s);

add_line("# alloc", s);
add_line(serialize(x->g_sw != nullptr), s);
if (x->g_sw) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->g_sw];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->g_sw[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# od_lw", s);

add_line("# alloc", s);
add_line(serialize(x->od_lw != nullptr), s);
if (x->od_lw) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->od_lw];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->od_lw[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# ssa_lw", s);

add_line("# alloc", s);
add_line(serialize(x->ssa_lw != nullptr), s);
if (x->ssa_lw) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->ssa_lw];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->ssa_lw[i]), s);
    }
}

}  // CONCLUDING IF
    std::string out = s.str();
    if (out.length() > 0) out.pop_back();
    return out;
}


std::string serialize(const gas_type* x) {
    std::stringstream s;
    
    std::string out = s.str();
    if (out.length() > 0) out.pop_back();
    return out;
}


std::string serialize(const config_type* x) {
    std::stringstream s;
    add_line("# i_band_from_reordered_g_lw", s);

add_line("# alloc", s);
add_line(serialize(x->i_band_from_reordered_g_lw != nullptr), s);
if (x->i_band_from_reordered_g_lw) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->i_band_from_reordered_g_lw];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->i_band_from_reordered_g_lw[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# i_band_from_reordered_g_sw", s);

add_line("# alloc", s);
add_line(serialize(x->i_band_from_reordered_g_sw != nullptr), s);
if (x->i_band_from_reordered_g_sw) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->i_band_from_reordered_g_sw];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->i_band_from_reordered_g_sw[i]), s);
    }
}

}  // CONCLUDING IF
    std::string out = s.str();
    if (out.length() > 0) out.pop_back();
    return out;
}


std::string serialize(const thermodynamics_type* x) {
    std::stringstream s;
    
    std::string out = s.str();
    if (out.length() > 0) out.pop_back();
    return out;
}

    
std::string config_injection(const aerosol_type& x) {
    std::stringstream out;
    out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_od_sw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_od_sw_d_0_s_19 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_od_sw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_od_sw_d_0_s_19 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_od_sw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_od_sw_d_1_s_20 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_od_sw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_od_sw_d_1_s_20 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_od_sw_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_od_sw_d_2_s_21 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_od_sw_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_od_sw_d_2_s_21 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_ssa_sw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_ssa_sw_d_0_s_22 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_ssa_sw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_ssa_sw_d_0_s_22 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_ssa_sw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_ssa_sw_d_1_s_23 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_ssa_sw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_ssa_sw_d_1_s_23 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_ssa_sw_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_ssa_sw_d_2_s_24 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_ssa_sw_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_ssa_sw_d_2_s_24 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_g_sw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_g_sw_d_0_s_25 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_g_sw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_g_sw_d_0_s_25 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_g_sw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_g_sw_d_1_s_26 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_g_sw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_g_sw_d_1_s_26 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_g_sw_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_g_sw_d_2_s_27 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_g_sw_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_g_sw_d_2_s_27 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_od_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_od_lw_d_0_s_28 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_od_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_od_lw_d_0_s_28 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_od_lw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_od_lw_d_1_s_29 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_od_lw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_od_lw_d_1_s_29 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_od_lw_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_od_lw_d_2_s_30 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_od_lw_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_od_lw_d_2_s_30 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_ssa_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_ssa_lw_d_0_s_31 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_ssa_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_ssa_lw_d_0_s_31 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_ssa_lw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_ssa_lw_d_1_s_32 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_ssa_lw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_ssa_lw_d_1_s_32 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SA_ssa_lw_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_ssa_lw_d_2_s_33 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"__f2dace_SOA_ssa_lw_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_ssa_lw_d_2_s_33 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"od_sw_a\", ";
out << "\"value\": \"" << (x.od_sw ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"ssa_sw_a\", ";
out << "\"value\": \"" << (x.ssa_sw ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"g_sw_a\", ";
out << "\"value\": \"" << (x.g_sw ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"od_lw_a\", ";
out << "\"value\": \"" << (x.od_lw ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_aerosol.aerosol_type\", ";
out << "\"component\": \"ssa_lw_a\", ";
out << "\"value\": \"" << (x.ssa_lw ? "true" : "false") << "\"}" << std::endl;
    return out.str();
}


std::string config_injection(const gas_type& x) {
    std::stringstream out;
    
    return out.str();
}


std::string config_injection(const config_type& x) {
    std::stringstream out;
    out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_config.config_type\", ";
out << "\"component\": \"__f2dace_SA_i_band_from_reordered_g_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_i_band_from_reordered_g_lw_d_0_s_2 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_config.config_type\", ";
out << "\"component\": \"__f2dace_SOA_i_band_from_reordered_g_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_2 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_config.config_type\", ";
out << "\"component\": \"__f2dace_SA_i_band_from_reordered_g_sw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_i_band_from_reordered_g_sw_d_0_s_3 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_config.config_type\", ";
out << "\"component\": \"__f2dace_SOA_i_band_from_reordered_g_sw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_3 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_config.config_type\", ";
out << "\"component\": \"i_band_from_reordered_g_lw_a\", ";
out << "\"value\": \"" << (x.i_band_from_reordered_g_lw ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_config.config_type\", ";
out << "\"component\": \"i_band_from_reordered_g_sw_a\", ";
out << "\"value\": \"" << (x.i_band_from_reordered_g_sw ? "true" : "false") << "\"}" << std::endl;
    return out.str();
}


std::string config_injection(const thermodynamics_type& x) {
    std::stringstream out;
    
    return out.str();
}

}  // namesepace serde

#endif // __DACE_SERDE__