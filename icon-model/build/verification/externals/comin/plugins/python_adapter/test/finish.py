# @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Please see the file LICENSE in the root of the source tree for this code.
# Where software is supplied by third parties, it is indicated in the
# headers of the routines.

import comin


@comin.register_callback(comin.EP_ATM_PHYSICS_BEFORE)
def phy():
    comin.finish("phy", "This is just a dummy error!")
