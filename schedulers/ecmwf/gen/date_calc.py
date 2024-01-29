#! /usr/bin/python

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

import optparse
import datetime
import calendar

#  -------------------------------------------------
#  DATE_CALC.PY
#  -------------------------------------------------
#
#  Python script doing some date calculations.
#
#  This script is part of the ICON SMS suite
#  Initial implementation: F. Prill, DWD (2012-05-07)
#
#  Corresponding author:
#    Florian Prill, DWD, mailto:florian.prill@dwd.de

# returns the last day of a month
def last_day_of_month(date):
    return (calendar.monthrange(date.year, date.month)[1])


def main():
    p = optparse.OptionParser()
    p.add_option('--date',     '-d', default="2012050100",  help="requested date YYYMMDDHH")
    p.add_option('--subtract', '-s', default="0",           help="days to subtract")
    p.add_option('--action',   '-a', default="printdate",   help="requested type of output: printdate/weekday/printfmt/endofmonth/fmt")
    p.add_option('--fmt',      '-f', default="",            help="requested date format (action=fmt)")
    #
    options, arguments = p.parse_args()

    # parse time string:
    if (options.date == "today"):
        date_object = datetime.datetime.now()
    else:
        date_object  = datetime.datetime.strptime(options.date, '%Y%m%d%H')
    # date calculation:
    date_object2 = date_object - datetime.timedelta(int(options.subtract))

    if (options.action == "printdate"):
        # print date (formatted output)
        print(date_object2.strftime("%Y%m%d%H")) 
    elif (options.action == "weekday"):
        # print only day of week
        print(date_object2.strftime("%w")) 
    elif (options.action == "printfmt"):
        # format as required by ICON:
        print(date_object2.strftime("%Y-%m-%dT%H:00:00Z"))
    elif (options.action == "fmt"):
        # format as required by ICON:
        print(date_object2.strftime(options.fmt))
    elif (options.action == "endofmonth"):
        # return "1" if end of month, "0" otherwise
        if (date_object2.day == last_day_of_month(date_object2)):
            print("1")
        else:
            print("0")
    else:
        print "Unknown command line option \"", options.action,"\"!"


if __name__ == '__main__':
    main()

