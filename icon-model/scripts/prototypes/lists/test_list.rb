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

$:.unshift File.join(File.dirname(__FILE__),".")
require 'list'
require 'minitest/autorun'
require 'pp'

OPERATORS = %w[mean min max div]
INTERVALS = %w[PT06H P01D P01M P01Y]
VARNAMES  = %w[u v temp zmc tend_t_up]

class MyListTest < Minitest::Test
  def setup
    @input = []
    OPERATORS.each {|operator| INTERVALS.each {|interval| VARNAMES.each {|varname|
          next if varname.size == 1 and interval == 'P01D'
          next if varname.size >  5 and interval == 'P01Y'
          next if varname.size == 4 and ['min','max'].include?(operator)
          @input << [operator,interval,varname]
    } } }
  end

  def intervalHasEvent?(interval)
    return interval == (ENV.has_key?('INT') ? ENV['INT'] : 'P01D')
  end

  # multiple adds of the same value, should lead to a singe entry only
  def check(listObj)
    listObj.clear
    listObj.add('a').add('a').add('a')
    assert_equal(1,listObj.size)
    listObj.add('b').add('b')
    assert_equal(2,listObj.size)
    pp listObj
  end
  def test_add_multiple
    check(MyHash.new)
    check(MyVector.new)
  end

  # Hash and Array based
  # OUTPUT = {
  # "intervalA" => [
  #    { "operatorA" => [vA,vB,...]},
  #    { "operatorB" => [vA,vC,...]},
  #    ],
  # "intervalB" => [
  #    { "operatorA" => [vA,vB,...]},
  #    { "operatorB" => [vD,vC,...]},
  #    ],
  # }
  def loop_hash(output)
    output.each {|interval,joblist|
      if intervalHasEvent?(interval) then
        puts "INTERVAL:#{interval} is active"
        joblist.each {|operator,varlist|
          varlist.each {|var|
            puts "Work on VARNAME:#{var} for OPERATOR:#{operator}"
          }
        }
      end
    }
  end
  def test_output_hash
    output = {}
    @input.each {|line|
      operator,interval,varname = line
      ((output[interval] ||= {})[operator] ||= []) << varname
    }
    loop_hash(output)
    pp output
  end

  # MyVector based
  $expectedOutput = MyVector.new([
    MyVector.new(["intervalA",
      MyVector.new(["operatorA", MyVector.new(['vA','vB'])]),
      MyVector.new(["operatorB", MyVector.new(['vA','vC'])])
    ]),
    MyVector.new(["intervalB",
      MyVector.new(["operatorA", MyVector.new(['vA','vB'])]),
      MyVector.new(["operatorB", MyVector.new(['vD','vC'])])
    ]),
  ])
  def loop_vector(output)
    puts "OUTPUT SIZE:#{output.size}"
    (0...output.size).each {|intervalIndex|
      intervalVector = output[intervalIndex]
      interval       = intervalVector[0]
      if intervalHasEvent?(interval) then
      else
        pp interval
      end
    }
  end
  def test_output_vector
    pp @input
    #pp $expectedOutput
    output         = MyVector.new
    intervalVector = MyVector.new
    operatorVector = MyVector.new
    varnameVector  = MyVector.new
    @input.each {|line|
      operator,interval,varname = line
      # check if interval is already  there
    }
    loop_vector($expectedOutput)
  end
end
