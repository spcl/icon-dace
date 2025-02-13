#ifndef __DACE_SERDE__
#define __DACE_SERDE__

#include <cassert>
#include <istream>
#include <iostream>
#include <sstream>
#include <optional>

#include "cloud_optics_fn_438.h"

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


void deserialize(cloud_type* x, std::istream& s) {
    bool yep;
    array_meta m;
    read_line(s, {"# mixing_ratio"});  // Should contain '# mixing_ratio'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_mixing_ratio_d_0_s_48 = m.size[0];
x->__f2dace_SA_mixing_ratio_d_1_s_49 = m.size[1];
x->__f2dace_SA_mixing_ratio_d_2_s_50 = m.size[2];
x->__f2dace_SOA_mixing_ratio_d_0_s_48 = m.lbound[0];
x->__f2dace_SOA_mixing_ratio_d_1_s_49 = m.lbound[1];
x->__f2dace_SOA_mixing_ratio_d_2_s_50 = m.lbound[2];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->mixing_ratio = new std::remove_pointer<decltype(x ->mixing_ratio)>::type[m.volume()];
ARRAY_META_DICT()[x->mixing_ratio] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->mixing_ratio[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# effective_radius"});  // Should contain '# effective_radius'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_effective_radius_d_0_s_51 = m.size[0];
x->__f2dace_SA_effective_radius_d_1_s_52 = m.size[1];
x->__f2dace_SA_effective_radius_d_2_s_53 = m.size[2];
x->__f2dace_SOA_effective_radius_d_0_s_51 = m.lbound[0];
x->__f2dace_SOA_effective_radius_d_1_s_52 = m.lbound[1];
x->__f2dace_SOA_effective_radius_d_2_s_53 = m.lbound[2];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->effective_radius = new std::remove_pointer<decltype(x ->effective_radius)>::type[m.volume()];
ARRAY_META_DICT()[x->effective_radius] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->effective_radius[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# fraction"});  // Should contain '# fraction'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_fraction_d_0_s_54 = m.size[0];
x->__f2dace_SA_fraction_d_1_s_55 = m.size[1];
x->__f2dace_SOA_fraction_d_0_s_54 = m.lbound[0];
x->__f2dace_SOA_fraction_d_1_s_55 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->fraction = new std::remove_pointer<decltype(x ->fraction)>::type[m.volume()];
ARRAY_META_DICT()[x->fraction] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->fraction[i]), s);
}

}  // CONCLUDING IF
}


void deserialize(cloud_optics_type* x, std::istream& s) {
    bool yep;
    array_meta m;
    read_line(s, {"# liq_coeff_lw"});  // Should contain '# liq_coeff_lw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_liq_coeff_lw_d_0_s_28 = m.size[0];
x->__f2dace_SA_liq_coeff_lw_d_1_s_29 = m.size[1];
x->__f2dace_SOA_liq_coeff_lw_d_0_s_28 = m.lbound[0];
x->__f2dace_SOA_liq_coeff_lw_d_1_s_29 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->liq_coeff_lw = new std::remove_pointer<decltype(x ->liq_coeff_lw)>::type[m.volume()];
ARRAY_META_DICT()[x->liq_coeff_lw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->liq_coeff_lw[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# liq_coeff_sw"});  // Should contain '# liq_coeff_sw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_liq_coeff_sw_d_0_s_30 = m.size[0];
x->__f2dace_SA_liq_coeff_sw_d_1_s_31 = m.size[1];
x->__f2dace_SOA_liq_coeff_sw_d_0_s_30 = m.lbound[0];
x->__f2dace_SOA_liq_coeff_sw_d_1_s_31 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->liq_coeff_sw = new std::remove_pointer<decltype(x ->liq_coeff_sw)>::type[m.volume()];
ARRAY_META_DICT()[x->liq_coeff_sw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->liq_coeff_sw[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# ice_coeff_lw"});  // Should contain '# ice_coeff_lw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_ice_coeff_lw_d_0_s_32 = m.size[0];
x->__f2dace_SA_ice_coeff_lw_d_1_s_33 = m.size[1];
x->__f2dace_SOA_ice_coeff_lw_d_0_s_32 = m.lbound[0];
x->__f2dace_SOA_ice_coeff_lw_d_1_s_33 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->ice_coeff_lw = new std::remove_pointer<decltype(x ->ice_coeff_lw)>::type[m.volume()];
ARRAY_META_DICT()[x->ice_coeff_lw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->ice_coeff_lw[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# ice_coeff_sw"});  // Should contain '# ice_coeff_sw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_ice_coeff_sw_d_0_s_34 = m.size[0];
x->__f2dace_SA_ice_coeff_sw_d_1_s_35 = m.size[1];
x->__f2dace_SOA_ice_coeff_sw_d_0_s_34 = m.lbound[0];
x->__f2dace_SOA_ice_coeff_sw_d_1_s_35 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->ice_coeff_sw = new std::remove_pointer<decltype(x ->ice_coeff_sw)>::type[m.volume()];
ARRAY_META_DICT()[x->ice_coeff_sw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->ice_coeff_sw[i]), s);
}

}  // CONCLUDING IF
}


void deserialize(config_type* x, std::istream& s) {
    bool yep;
    array_meta m;
    read_line(s, {"# cloud_optics"});  // Should contain '# cloud_optics'

x ->cloud_optics = new std::remove_pointer<decltype(x ->cloud_optics)>::type;
deserialize(x->cloud_optics, s);

}


void deserialize(thermodynamics_type* x, std::istream& s) {
    bool yep;
    array_meta m;
    read_line(s, {"# pressure_hl"});  // Should contain '# pressure_hl'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_pressure_hl_d_0_s_38 = m.size[0];
x->__f2dace_SA_pressure_hl_d_1_s_39 = m.size[1];
x->__f2dace_SOA_pressure_hl_d_0_s_38 = m.lbound[0];
x->__f2dace_SOA_pressure_hl_d_1_s_39 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->pressure_hl = new std::remove_pointer<decltype(x ->pressure_hl)>::type[m.volume()];
ARRAY_META_DICT()[x->pressure_hl] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->pressure_hl[i]), s);
}

}  // CONCLUDING IF
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


std::string serialize(const cloud_type* x) {
    std::stringstream s;
    add_line("# mixing_ratio", s);

add_line("# alloc", s);
add_line(serialize(x->mixing_ratio != nullptr), s);
if (x->mixing_ratio) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->mixing_ratio];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->mixing_ratio[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# effective_radius", s);

add_line("# alloc", s);
add_line(serialize(x->effective_radius != nullptr), s);
if (x->effective_radius) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->effective_radius];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->effective_radius[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# fraction", s);

add_line("# alloc", s);
add_line(serialize(x->fraction != nullptr), s);
if (x->fraction) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->fraction];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->fraction[i]), s);
    }
}

}  // CONCLUDING IF
    std::string out = s.str();
    if (out.length() > 0) out.pop_back();
    return out;
}


std::string serialize(const cloud_optics_type* x) {
    std::stringstream s;
    add_line("# liq_coeff_lw", s);

add_line("# alloc", s);
add_line(serialize(x->liq_coeff_lw != nullptr), s);
if (x->liq_coeff_lw) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->liq_coeff_lw];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->liq_coeff_lw[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# liq_coeff_sw", s);

add_line("# alloc", s);
add_line(serialize(x->liq_coeff_sw != nullptr), s);
if (x->liq_coeff_sw) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->liq_coeff_sw];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->liq_coeff_sw[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# ice_coeff_lw", s);

add_line("# alloc", s);
add_line(serialize(x->ice_coeff_lw != nullptr), s);
if (x->ice_coeff_lw) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->ice_coeff_lw];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->ice_coeff_lw[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# ice_coeff_sw", s);

add_line("# alloc", s);
add_line(serialize(x->ice_coeff_sw != nullptr), s);
if (x->ice_coeff_sw) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->ice_coeff_sw];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->ice_coeff_sw[i]), s);
    }
}

}  // CONCLUDING IF
    std::string out = s.str();
    if (out.length() > 0) out.pop_back();
    return out;
}


std::string serialize(const config_type* x) {
    std::stringstream s;
    add_line("# cloud_optics", s);
add_line(serialize(x->cloud_optics), s);
    std::string out = s.str();
    if (out.length() > 0) out.pop_back();
    return out;
}


std::string serialize(const thermodynamics_type* x) {
    std::stringstream s;
    add_line("# pressure_hl", s);

add_line("# alloc", s);
add_line(serialize(x->pressure_hl != nullptr), s);
if (x->pressure_hl) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->pressure_hl];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->pressure_hl[i]), s);
    }
}

}  // CONCLUDING IF
    std::string out = s.str();
    if (out.length() > 0) out.pop_back();
    return out;
}

    
std::string config_injection(const cloud_type& x) {
    std::stringstream out;
    out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_mixing_ratio_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_mixing_ratio_d_0_s_48 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_mixing_ratio_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_mixing_ratio_d_0_s_48 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_mixing_ratio_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_mixing_ratio_d_1_s_49 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_mixing_ratio_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_mixing_ratio_d_1_s_49 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_mixing_ratio_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_mixing_ratio_d_2_s_50 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_mixing_ratio_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_mixing_ratio_d_2_s_50 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_effective_radius_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_effective_radius_d_0_s_51 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_effective_radius_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_effective_radius_d_0_s_51 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_effective_radius_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_effective_radius_d_1_s_52 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_effective_radius_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_effective_radius_d_1_s_52 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_effective_radius_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_effective_radius_d_2_s_53 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_effective_radius_d_2_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_effective_radius_d_2_s_53 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_fraction_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_fraction_d_0_s_54 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_fraction_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_fraction_d_0_s_54 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_fraction_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_fraction_d_1_s_55 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_fraction_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_fraction_d_1_s_55 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"mixing_ratio_a\", ";
out << "\"value\": \"" << (x.mixing_ratio ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"effective_radius_a\", ";
out << "\"value\": \"" << (x.effective_radius ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"fraction_a\", ";
out << "\"value\": \"" << (x.fraction ? "true" : "false") << "\"}" << std::endl;
    return out.str();
}


std::string config_injection(const cloud_optics_type& x) {
    std::stringstream out;
    out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SA_liq_coeff_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_liq_coeff_lw_d_0_s_28 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SOA_liq_coeff_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_liq_coeff_lw_d_0_s_28 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SA_liq_coeff_lw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_liq_coeff_lw_d_1_s_29 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SOA_liq_coeff_lw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_liq_coeff_lw_d_1_s_29 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SA_liq_coeff_sw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_liq_coeff_sw_d_0_s_30 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SOA_liq_coeff_sw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_liq_coeff_sw_d_0_s_30 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SA_liq_coeff_sw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_liq_coeff_sw_d_1_s_31 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SOA_liq_coeff_sw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_liq_coeff_sw_d_1_s_31 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SA_ice_coeff_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_ice_coeff_lw_d_0_s_32 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SOA_ice_coeff_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_ice_coeff_lw_d_0_s_32 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SA_ice_coeff_lw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_ice_coeff_lw_d_1_s_33 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SOA_ice_coeff_lw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_ice_coeff_lw_d_1_s_33 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SA_ice_coeff_sw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_ice_coeff_sw_d_0_s_34 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SOA_ice_coeff_sw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_ice_coeff_sw_d_0_s_34 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SA_ice_coeff_sw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_ice_coeff_sw_d_1_s_35 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"__f2dace_SOA_ice_coeff_sw_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_ice_coeff_sw_d_1_s_35 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"liq_coeff_lw_a\", ";
out << "\"value\": \"" << (x.liq_coeff_lw ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"liq_coeff_sw_a\", ";
out << "\"value\": \"" << (x.liq_coeff_sw ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"ice_coeff_lw_a\", ";
out << "\"value\": \"" << (x.ice_coeff_lw ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud_optics_data.cloud_optics_type\", ";
out << "\"component\": \"ice_coeff_sw_a\", ";
out << "\"value\": \"" << (x.ice_coeff_sw ? "true" : "false") << "\"}" << std::endl;
    return out.str();
}


std::string config_injection(const config_type& x) {
    std::stringstream out;
    
    return out.str();
}


std::string config_injection(const thermodynamics_type& x) {
    std::stringstream out;
    out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_thermodynamics.thermodynamics_type\", ";
out << "\"component\": \"__f2dace_SA_pressure_hl_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_pressure_hl_d_0_s_38 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_thermodynamics.thermodynamics_type\", ";
out << "\"component\": \"__f2dace_SOA_pressure_hl_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_pressure_hl_d_0_s_38 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_thermodynamics.thermodynamics_type\", ";
out << "\"component\": \"__f2dace_SA_pressure_hl_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_pressure_hl_d_1_s_39 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_thermodynamics.thermodynamics_type\", ";
out << "\"component\": \"__f2dace_SOA_pressure_hl_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_pressure_hl_d_1_s_39 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_thermodynamics.thermodynamics_type\", ";
out << "\"component\": \"pressure_hl_a\", ";
out << "\"value\": \"" << (x.pressure_hl ? "true" : "false") << "\"}" << std::endl;
    return out.str();
}

}  // namesepace serde

#endif // __DACE_SERDE__