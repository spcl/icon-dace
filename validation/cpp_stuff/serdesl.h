#ifndef __DACE_SERDE__
#define __DACE_SERDE__

#include <cassert>
#include <istream>
#include <iostream>
#include <sstream>
#include <optional>

#include "solver_mcica_lw.h"

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
    read_line(s, {"# fraction"});  // Should contain '# fraction'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_fraction_d_0_s_42 = m.size[0];
x->__f2dace_SA_fraction_d_1_s_43 = m.size[1];
x->__f2dace_SOA_fraction_d_0_s_42 = m.lbound[0];
x->__f2dace_SOA_fraction_d_1_s_43 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->fraction = new std::remove_pointer<decltype(x ->fraction)>::type[m.volume()];
ARRAY_META_DICT()[x->fraction] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->fraction[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# fractional_std"});  // Should contain '# fractional_std'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_fractional_std_d_0_s_44 = m.size[0];
x->__f2dace_SA_fractional_std_d_1_s_45 = m.size[1];
x->__f2dace_SOA_fractional_std_d_0_s_44 = m.lbound[0];
x->__f2dace_SOA_fractional_std_d_1_s_45 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->fractional_std = new std::remove_pointer<decltype(x ->fractional_std)>::type[m.volume()];
ARRAY_META_DICT()[x->fractional_std] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->fractional_std[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# overlap_param"});  // Should contain '# overlap_param'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_overlap_param_d_0_s_46 = m.size[0];
x->__f2dace_SA_overlap_param_d_1_s_47 = m.size[1];
x->__f2dace_SOA_overlap_param_d_0_s_46 = m.lbound[0];
x->__f2dace_SOA_overlap_param_d_1_s_47 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->overlap_param = new std::remove_pointer<decltype(x ->overlap_param)>::type[m.volume()];
ARRAY_META_DICT()[x->overlap_param] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->overlap_param[i]), s);
}

}  // CONCLUDING IF
}


void deserialize(pdf_sampler_type* x, std::istream& s) {
    bool yep;
    array_meta m;
    read_line(s, {"# ncdf"});  // Should contain '# ncdf'

deserialize(&(x->ncdf), s);

read_line(s, {"# nfsd"});  // Should contain '# nfsd'

deserialize(&(x->nfsd), s);

read_line(s, {"# fsd1"});  // Should contain '# fsd1'

deserialize(&(x->fsd1), s);

read_line(s, {"# inv_fsd_interval"});  // Should contain '# inv_fsd_interval'

deserialize(&(x->inv_fsd_interval), s);

read_line(s, {"# val"});  // Should contain '# val'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_val_d_0_s_4 = m.size[0];
x->__f2dace_SA_val_d_1_s_5 = m.size[1];
x->__f2dace_SOA_val_d_0_s_4 = m.lbound[0];
x->__f2dace_SOA_val_d_1_s_5 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->val = new std::remove_pointer<decltype(x ->val)>::type[m.volume()];
ARRAY_META_DICT()[x->val] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->val[i]), s);
}

}  // CONCLUDING IF
}


void deserialize(config_type* x, std::istream& s) {
    bool yep;
    array_meta m;
    read_line(s, {"# i_band_from_reordered_g_lw"});  // Should contain '# i_band_from_reordered_g_lw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_i_band_from_reordered_g_lw_d_0_s_7 = m.size[0];
x->__f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_7 = m.lbound[0];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->i_band_from_reordered_g_lw = new std::remove_pointer<decltype(x ->i_band_from_reordered_g_lw)>::type[m.volume()];
ARRAY_META_DICT()[x->i_band_from_reordered_g_lw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->i_band_from_reordered_g_lw[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# pdf_sampler"});  // Should contain '# pdf_sampler'

x ->pdf_sampler = new std::remove_pointer<decltype(x ->pdf_sampler)>::type;
deserialize(x->pdf_sampler, s);

}


void deserialize(flux_type* x, std::istream& s) {
    bool yep;
    array_meta m;
    read_line(s, {"# lw_up"});  // Should contain '# lw_up'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_lw_up_d_0_s_21 = m.size[0];
x->__f2dace_SA_lw_up_d_1_s_22 = m.size[1];
x->__f2dace_SOA_lw_up_d_0_s_21 = m.lbound[0];
x->__f2dace_SOA_lw_up_d_1_s_22 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->lw_up = new std::remove_pointer<decltype(x ->lw_up)>::type[m.volume()];
ARRAY_META_DICT()[x->lw_up] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->lw_up[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# lw_dn"});  // Should contain '# lw_dn'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_lw_dn_d_0_s_23 = m.size[0];
x->__f2dace_SA_lw_dn_d_1_s_24 = m.size[1];
x->__f2dace_SOA_lw_dn_d_0_s_23 = m.lbound[0];
x->__f2dace_SOA_lw_dn_d_1_s_24 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->lw_dn = new std::remove_pointer<decltype(x ->lw_dn)>::type[m.volume()];
ARRAY_META_DICT()[x->lw_dn] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->lw_dn[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# lw_up_clear"});  // Should contain '# lw_up_clear'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_lw_up_clear_d_0_s_25 = m.size[0];
x->__f2dace_SA_lw_up_clear_d_1_s_26 = m.size[1];
x->__f2dace_SOA_lw_up_clear_d_0_s_25 = m.lbound[0];
x->__f2dace_SOA_lw_up_clear_d_1_s_26 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->lw_up_clear = new std::remove_pointer<decltype(x ->lw_up_clear)>::type[m.volume()];
ARRAY_META_DICT()[x->lw_up_clear] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->lw_up_clear[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# lw_dn_clear"});  // Should contain '# lw_dn_clear'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_lw_dn_clear_d_0_s_27 = m.size[0];
x->__f2dace_SA_lw_dn_clear_d_1_s_28 = m.size[1];
x->__f2dace_SOA_lw_dn_clear_d_0_s_27 = m.lbound[0];
x->__f2dace_SOA_lw_dn_clear_d_1_s_28 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->lw_dn_clear = new std::remove_pointer<decltype(x ->lw_dn_clear)>::type[m.volume()];
ARRAY_META_DICT()[x->lw_dn_clear] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->lw_dn_clear[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# lw_dn_surf_g"});  // Should contain '# lw_dn_surf_g'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_lw_dn_surf_g_d_0_s_29 = m.size[0];
x->__f2dace_SA_lw_dn_surf_g_d_1_s_30 = m.size[1];
x->__f2dace_SOA_lw_dn_surf_g_d_0_s_29 = m.lbound[0];
x->__f2dace_SOA_lw_dn_surf_g_d_1_s_30 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->lw_dn_surf_g = new std::remove_pointer<decltype(x ->lw_dn_surf_g)>::type[m.volume()];
ARRAY_META_DICT()[x->lw_dn_surf_g] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->lw_dn_surf_g[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# lw_dn_surf_clear_g"});  // Should contain '# lw_dn_surf_clear_g'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_lw_dn_surf_clear_g_d_0_s_31 = m.size[0];
x->__f2dace_SA_lw_dn_surf_clear_g_d_1_s_32 = m.size[1];
x->__f2dace_SOA_lw_dn_surf_clear_g_d_0_s_31 = m.lbound[0];
x->__f2dace_SOA_lw_dn_surf_clear_g_d_1_s_32 = m.lbound[1];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->lw_dn_surf_clear_g = new std::remove_pointer<decltype(x ->lw_dn_surf_clear_g)>::type[m.volume()];
ARRAY_META_DICT()[x->lw_dn_surf_clear_g] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->lw_dn_surf_clear_g[i]), s);
}

}  // CONCLUDING IF
read_line(s, {"# cloud_cover_lw"});  // Should contain '# cloud_cover_lw'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_cloud_cover_lw_d_0_s_33 = m.size[0];
x->__f2dace_SOA_cloud_cover_lw_d_0_s_33 = m.lbound[0];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->cloud_cover_lw = new std::remove_pointer<decltype(x ->cloud_cover_lw)>::type[m.volume()];
ARRAY_META_DICT()[x->cloud_cover_lw] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->cloud_cover_lw[i]), s);
}

}  // CONCLUDING IF
}


void deserialize(randomnumberstream* x, std::istream& s) {
    bool yep;
    array_meta m;
    read_line(s, {"# iused"});  // Should contain '# iused'

deserialize(&(x->iused), s);

read_line(s, {"# inittest"});  // Should contain '# inittest'

deserialize(&(x->inittest), s);

read_line(s, {"# ix"});  // Should contain '# ix'

m = read_array_meta(s);


read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->ix = new std::remove_pointer<decltype(x ->ix)>::type[m.volume()];
ARRAY_META_DICT()[x->ix] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->ix[i]), s);
}

read_line(s, {"# zrm"});  // Should contain '# zrm'

deserialize(&(x->zrm), s);

}


void deserialize(single_level_type* x, std::istream& s) {
    bool yep;
    array_meta m;
    read_line(s, {"# iseed"});  // Should contain '# iseed'

read_line(s, {"# alloc"});  // Should contain '# alloc'
deserialize(&yep, s);
if (yep) {  // BEGINING IF


m = read_array_meta(s);
x->__f2dace_SA_iseed_d_0_s_35 = m.size[0];
x->__f2dace_SOA_iseed_d_0_s_35 = m.lbound[0];
read_line(s, {"# entries"});  // Should contain '# entries'
// We only need to allocate a volume of contiguous memory, and let DaCe interpret (assuming it follows the same protocol 
// as us).
x ->iseed = new std::remove_pointer<decltype(x ->iseed)>::type[m.volume()];
ARRAY_META_DICT()[x->iseed] = m;
for (int i=0; i<m.volume(); ++i) {
  deserialize(&(x->iseed[i]), s);
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
add_line("# fractional_std", s);

add_line("# alloc", s);
add_line(serialize(x->fractional_std != nullptr), s);
if (x->fractional_std) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->fractional_std];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->fractional_std[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# overlap_param", s);

add_line("# alloc", s);
add_line(serialize(x->overlap_param != nullptr), s);
if (x->overlap_param) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->overlap_param];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->overlap_param[i]), s);
    }
}

}  // CONCLUDING IF
    std::string out = s.str();
    if (out.length() > 0) out.pop_back();
    return out;
}


std::string serialize(const pdf_sampler_type* x) {
    std::stringstream s;
    add_line("# ncdf", s);
add_line(serialize(x->ncdf), s);
add_line("# nfsd", s);
add_line(serialize(x->nfsd), s);
add_line("# fsd1", s);
add_line(serialize(x->fsd1), s);
add_line("# inv_fsd_interval", s);
add_line(serialize(x->inv_fsd_interval), s);
add_line("# val", s);

add_line("# alloc", s);
add_line(serialize(x->val != nullptr), s);
if (x->val) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->val];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->val[i]), s);
    }
}

}  // CONCLUDING IF
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
add_line("# pdf_sampler", s);
add_line(serialize(x->pdf_sampler), s);
    std::string out = s.str();
    if (out.length() > 0) out.pop_back();
    return out;
}


std::string serialize(const flux_type* x) {
    std::stringstream s;
    add_line("# lw_up", s);

add_line("# alloc", s);
add_line(serialize(x->lw_up != nullptr), s);
if (x->lw_up) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->lw_up];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->lw_up[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# lw_dn", s);

add_line("# alloc", s);
add_line(serialize(x->lw_dn != nullptr), s);
if (x->lw_dn) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->lw_dn];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->lw_dn[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# lw_up_clear", s);

add_line("# alloc", s);
add_line(serialize(x->lw_up_clear != nullptr), s);
if (x->lw_up_clear) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->lw_up_clear];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->lw_up_clear[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# lw_dn_clear", s);

add_line("# alloc", s);
add_line(serialize(x->lw_dn_clear != nullptr), s);
if (x->lw_dn_clear) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->lw_dn_clear];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->lw_dn_clear[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# lw_dn_surf_g", s);

add_line("# alloc", s);
add_line(serialize(x->lw_dn_surf_g != nullptr), s);
if (x->lw_dn_surf_g) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->lw_dn_surf_g];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->lw_dn_surf_g[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# lw_dn_surf_clear_g", s);

add_line("# alloc", s);
add_line(serialize(x->lw_dn_surf_clear_g != nullptr), s);
if (x->lw_dn_surf_clear_g) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->lw_dn_surf_clear_g];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->lw_dn_surf_clear_g[i]), s);
    }
}

}  // CONCLUDING IF
add_line("# cloud_cover_lw", s);

add_line("# alloc", s);
add_line(serialize(x->cloud_cover_lw != nullptr), s);
if (x->cloud_cover_lw) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->cloud_cover_lw];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->cloud_cover_lw[i]), s);
    }
}

}  // CONCLUDING IF
    std::string out = s.str();
    if (out.length() > 0) out.pop_back();
    return out;
}


std::string serialize(const randomnumberstream* x) {
    std::stringstream s;
    add_line("# iused", s);
add_line(serialize(x->iused), s);
add_line("# inittest", s);
add_line(serialize(x->inittest), s);
add_line("# ix", s);

{
    const array_meta& m = ARRAY_META_DICT()[x->ix];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->ix[i]), s);
    }
}

add_line("# zrm", s);
add_line(serialize(x->zrm), s);
    std::string out = s.str();
    if (out.length() > 0) out.pop_back();
    return out;
}


std::string serialize(const single_level_type* x) {
    std::stringstream s;
    add_line("# iseed", s);

add_line("# alloc", s);
add_line(serialize(x->iseed != nullptr), s);
if (x->iseed) {    // BEGINING IF


{
    const array_meta& m = ARRAY_META_DICT()[x->iseed];
    add_line("# rank", s);
    add_line(m.rank, s);
    add_line("# size", s);
    for (auto i : m.size) add_line(i, s);
    add_line("# lbound", s);
    for (auto i : m.lbound) add_line(i, s);
    add_line("# entries", s);
    for (int i=0; i<m.volume(); ++i) {
        add_line(serialize(x->iseed[i]), s);
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
out << "\"component\": \"__f2dace_SA_fraction_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_fraction_d_0_s_42 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_fraction_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_fraction_d_0_s_42 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_fraction_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_fraction_d_1_s_43 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_fraction_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_fraction_d_1_s_43 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_fractional_std_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_fractional_std_d_0_s_44 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_fractional_std_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_fractional_std_d_0_s_44 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_fractional_std_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_fractional_std_d_1_s_45 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_fractional_std_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_fractional_std_d_1_s_45 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_overlap_param_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_overlap_param_d_0_s_46 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_overlap_param_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_overlap_param_d_0_s_46 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SA_overlap_param_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_overlap_param_d_1_s_47 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"__f2dace_SOA_overlap_param_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_overlap_param_d_1_s_47 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"fraction_a\", ";
out << "\"value\": \"" << (x.fraction ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"fractional_std_a\", ";
out << "\"value\": \"" << (x.fractional_std ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_cloud.cloud_type\", ";
out << "\"component\": \"overlap_param_a\", ";
out << "\"value\": \"" << (x.overlap_param ? "true" : "false") << "\"}" << std::endl;
    return out.str();
}


std::string config_injection(const pdf_sampler_type& x) {
    std::stringstream out;
    out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_pdf_sampler.pdf_sampler_type\", ";
out << "\"component\": \"ncdf\", ";
out << "\"value\": \"" << x.ncdf << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_pdf_sampler.pdf_sampler_type\", ";
out << "\"component\": \"nfsd\", ";
out << "\"value\": \"" << x.nfsd << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_pdf_sampler.pdf_sampler_type\", ";
out << "\"component\": \"fsd1\", ";
out << "\"value\": \"" << x.fsd1 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_pdf_sampler.pdf_sampler_type\", ";
out << "\"component\": \"inv_fsd_interval\", ";
out << "\"value\": \"" << x.inv_fsd_interval << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_pdf_sampler.pdf_sampler_type\", ";
out << "\"component\": \"__f2dace_SA_val_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_val_d_0_s_4 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_pdf_sampler.pdf_sampler_type\", ";
out << "\"component\": \"__f2dace_SOA_val_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_val_d_0_s_4 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_pdf_sampler.pdf_sampler_type\", ";
out << "\"component\": \"__f2dace_SA_val_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_val_d_1_s_5 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_pdf_sampler.pdf_sampler_type\", ";
out << "\"component\": \"__f2dace_SOA_val_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_val_d_1_s_5 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_pdf_sampler.pdf_sampler_type\", ";
out << "\"component\": \"val_a\", ";
out << "\"value\": \"" << (x.val ? "true" : "false") << "\"}" << std::endl;
    return out.str();
}


std::string config_injection(const config_type& x) {
    std::stringstream out;
    out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_config.config_type\", ";
out << "\"component\": \"__f2dace_SA_i_band_from_reordered_g_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_i_band_from_reordered_g_lw_d_0_s_7 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_config.config_type\", ";
out << "\"component\": \"__f2dace_SOA_i_band_from_reordered_g_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_7 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_config.config_type\", ";
out << "\"component\": \"i_band_from_reordered_g_lw_a\", ";
out << "\"value\": \"" << (x.i_band_from_reordered_g_lw ? "true" : "false") << "\"}" << std::endl;
    return out.str();
}


std::string config_injection(const flux_type& x) {
    std::stringstream out;
    out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_lw_up_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_lw_up_d_0_s_21 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_lw_up_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_lw_up_d_0_s_21 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_lw_up_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_lw_up_d_1_s_22 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_lw_up_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_lw_up_d_1_s_22 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_lw_dn_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_lw_dn_d_0_s_23 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_lw_dn_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_lw_dn_d_0_s_23 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_lw_dn_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_lw_dn_d_1_s_24 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_lw_dn_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_lw_dn_d_1_s_24 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_lw_up_clear_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_lw_up_clear_d_0_s_25 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_lw_up_clear_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_lw_up_clear_d_0_s_25 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_lw_up_clear_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_lw_up_clear_d_1_s_26 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_lw_up_clear_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_lw_up_clear_d_1_s_26 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_lw_dn_clear_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_lw_dn_clear_d_0_s_27 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_lw_dn_clear_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_lw_dn_clear_d_0_s_27 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_lw_dn_clear_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_lw_dn_clear_d_1_s_28 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_lw_dn_clear_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_lw_dn_clear_d_1_s_28 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_lw_dn_surf_g_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_lw_dn_surf_g_d_0_s_29 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_lw_dn_surf_g_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_lw_dn_surf_g_d_0_s_29 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_lw_dn_surf_g_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_lw_dn_surf_g_d_1_s_30 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_lw_dn_surf_g_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_lw_dn_surf_g_d_1_s_30 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_lw_dn_surf_clear_g_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_lw_dn_surf_clear_g_d_0_s_31 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_lw_dn_surf_clear_g_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_lw_dn_surf_clear_g_d_0_s_31 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_lw_dn_surf_clear_g_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_lw_dn_surf_clear_g_d_1_s_32 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_lw_dn_surf_clear_g_d_1_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_lw_dn_surf_clear_g_d_1_s_32 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SA_cloud_cover_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_cloud_cover_lw_d_0_s_33 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"__f2dace_SOA_cloud_cover_lw_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_cloud_cover_lw_d_0_s_33 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"lw_up_a\", ";
out << "\"value\": \"" << (x.lw_up ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"lw_dn_a\", ";
out << "\"value\": \"" << (x.lw_dn ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"lw_up_clear_a\", ";
out << "\"value\": \"" << (x.lw_up_clear ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"lw_dn_clear_a\", ";
out << "\"value\": \"" << (x.lw_dn_clear ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"lw_dn_surf_g_a\", ";
out << "\"value\": \"" << (x.lw_dn_surf_g ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"lw_dn_surf_clear_g_a\", ";
out << "\"value\": \"" << (x.lw_dn_surf_clear_g ? "true" : "false") << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_flux.flux_type\", ";
out << "\"component\": \"cloud_cover_lw_a\", ";
out << "\"value\": \"" << (x.cloud_cover_lw ? "true" : "false") << "\"}" << std::endl;
    return out.str();
}


std::string config_injection(const randomnumberstream& x) {
    std::stringstream out;
    out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"random_numbers_mix.randomnumberstream\", ";
out << "\"component\": \"iused\", ";
out << "\"value\": \"" << x.iused << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"random_numbers_mix.randomnumberstream\", ";
out << "\"component\": \"inittest\", ";
out << "\"value\": \"" << x.inittest << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"random_numbers_mix.randomnumberstream\", ";
out << "\"component\": \"zrm\", ";
out << "\"value\": \"" << x.zrm << "\"}" << std::endl;
    return out.str();
}


std::string config_injection(const single_level_type& x) {
    std::stringstream out;
    out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_single_level.single_level_type\", ";
out << "\"component\": \"__f2dace_SA_iseed_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SA_iseed_d_0_s_35 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_single_level.single_level_type\", ";
out << "\"component\": \"__f2dace_SOA_iseed_d_0_s\", ";
out << "\"value\": \"" << x.__f2dace_SOA_iseed_d_0_s_35 << "\"}" << std::endl;
out << "{";
out << "\"type\": \"ConstTypeInjection\", ";
out << "\"scope\": null, ";
out << "\"root\": \"radiation_single_level.single_level_type\", ";
out << "\"component\": \"iseed_a\", ";
out << "\"value\": \"" << (x.iseed ? "true" : "false") << "\"}" << std::endl;
    return out.str();
}

}  // namesepace serde

#endif // __DACE_SERDE__