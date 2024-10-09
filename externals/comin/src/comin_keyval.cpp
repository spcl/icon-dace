// @file comin_keyval.cpp
// @brief Key value storage backend based on std::unordered_map and std::variant
//
// @authors 08/2024 :: ICON Community Interface  <comin@icon-model.org>
//
// SPDX-License-Identifier: BSD-3-Clause
//
// See LICENSES for license information.
// Where software is supplied by third parties, it is indicated in the
// headers of the routines.

#include <iostream>
#include <variant>
#include <string>
#include <unordered_map>

namespace {
  using VariantType = std::variant<int, double, std::string, bool>;
  using MapType = std::unordered_map<std::string, VariantType>;
}

extern "C" {
  void comin_keyval_set_int_c(const char* ckey, int val, MapType* map){
    (*map)[ckey] = val;
  };

  void comin_keyval_get_int_c(const char* ckey, int* val, MapType* map){
    *val = std::get<int>(map->at(ckey));
  };

  void comin_keyval_set_double_c(const char* ckey, double val, MapType* map){
    (*map)[ckey] = val;
  };

  void comin_keyval_get_double_c(const char* ckey, double* val, MapType* map){
    *val = std::get<double>(map->at(ckey));
  };

  void comin_keyval_set_char_c(const char* ckey, char* val, MapType* map){
    (*map)[ckey] = std::string(val);
  };

  void comin_keyval_get_char_c(const char* ckey, const char** val, MapType* map){
    *val = std::get<std::string>(map->at(ckey)).data();
  };

  void comin_keyval_set_bool_c(const char* ckey, bool val, MapType* map){
    (*map)[ckey] = val;
  };

  void comin_keyval_get_bool_c(const char* ckey, bool *val, MapType* map){
    *val = std::get<bool>(map->at(ckey));
  };

  void comin_keyval_create_c(MapType ** map)
  {
    *map = new MapType();
  };

  void comin_keyval_delete_c(MapType * map){
    delete map;
  };

  void comin_keyval_query_c(const char* ckey, int* idx, MapType * map){
    auto it = map->find(ckey);
    if (it == map->end()) *idx = -1;
    else *idx = it->second.index();
  };

  void comin_keyval_iterator_begin_c(MapType * map, MapType::iterator** it){
    *it = new MapType::iterator(map->begin());
  }

  void comin_keyval_iterator_end_c(MapType * map, MapType::iterator** it){
    *it = new MapType::iterator(map->end());
  }

  const char* comin_keyval_iterator_get_key_c(MapType::iterator* it){
    return ((*it)->first).c_str();
  }

  bool comin_keyval_iterator_compare_c(MapType::iterator* it1, MapType::iterator* it2){
    return(*it1 == *it2);
  }

  void comin_keyval_iterator_next_c(MapType::iterator* it){
    (*it)++;
  }

  void comin_keyval_iterator_delete_c(MapType::iterator* it){
    delete it;
  }

}
