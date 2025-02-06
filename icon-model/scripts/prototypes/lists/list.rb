# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

# Author: ralf.mueller

class MyVector < Array
  def add(value)
    self << value unless self.include?(value)
    self
  end
end
class MyHash < Hash
  def add(value)
    self[value] = {} unless self.keys.include?(value)
    self
  end
end
