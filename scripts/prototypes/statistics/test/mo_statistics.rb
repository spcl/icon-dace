# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

$:.unshift File.join(File.dirname(__FILE__),"..",".")
require "test/unit"
require "mo_statistics"
require "pp"

class TestCodeParser < Test::Unit::TestCase
  def test_create_single_stat
    # check plain scalar contructor
    (1..4).each {|i|
      ds = DataStatistics.new
      assert_equal([],ds.accumulation)
    }
  end
  def test_stat_with_target
    # check contructor with target
    target = []
    ds = DataStatistics.new(target: target)
    assert_equal(target, ds.accumulation)

    ds.add(NArray[[1.0,1.0],[1.0,1.0]])
    assert_equal(NArray[[1.0,1.0],[1.0,1.0]],ds.accumulation[0])

    ds.add(NArray[[1.0,1.0],[1.0,1.0]])
    assert_equal(NArray[[2.0,2.0],[2.0,2.0]],ds.accumulation[0])

    assert_equal(target,ds.accumulation) 
  end
  def test_stat_with_source
    # check contructor with source
    source = [NArray[[1.0,2.0],[3.0,4.0]]]
    ds = DataStatistics.new(source: source)

    ds.add
    assert_equal(source, ds.accumulation)
    ds.add
    assert_equal(2*source[0], ds.accumulation[0])
    ds.add
    ds.add
    assert_equal(4*source[0], ds.accumulation[0])

    source[0] = -NArray[[1.0,2.0],[3.0,4.0]]
    ds.add
    assert_equal(-3*source[0], ds.accumulation[0])
  end
  def test_stat_with_source_and_target
    source = [NArray.float(2,2)+1]
    target = [NArray.float(2,2)+10]
    ds = DataStatistics.new(source: source, target: target)

    ds.add
    assert_equal([NArray.float(2,2)+11],ds.accumulation)

    source[0] = NArray.float(2,2)+2
    ds.add
    assert_equal([NArray.float(2,2)+13],ds.accumulation)
  end
  def test_init
    test_create_single_stat
    test_stat_with_target
    test_stat_with_source
    test_stat_with_source_and_target
  end
  def test_addTo_single_stat
    addMe = NArray[1.0,2.0,3.0]
    ds = DataStatistics.new
    ds.add(addMe)
    assert_equal(1,ds.numberOfAccumulations)
    ds.add(addMe)
    assert_equal(2,ds.numberOfAccumulations)
    10.times { ds.add(addMe) }
    assert_equal(12,ds.numberOfAccumulations)
    assert_equal(12.0,ds.accumulation.first[0])
    assert_equal(24.0,ds.accumulation.first[1])
    assert_equal(36.0,ds.accumulation.first[2])

    assert_equal(1.0,ds.min.uniq[0])
    assert_equal(2.0,ds.mean.uniq[0])
    assert_equal(3.0,ds.max.uniq[0])
    assert_equal(2.0,ds.median.uniq[0])
  end
  def test_add_usage
    ds = DataStatistics.new
    assert_raise ArgumentError do; ds.add; end
    ds = DataStatistics.new(source: [1])
    assert_raise ArgumentError do; ds.add(2); end
  end

  def test_stat_collections
    oceTemp       = [NArray[24.0]]
    ds_oceTemp    = DataStatistics.new(source: oceTemp)

    oceSal        = [NArray[35.0]]
    ds_oceSal     = DataStatistics.new(source: oceSal)

    oceDens       = [NArray[1029.0]]
    ds_oceDens    = DataStatistics.new(source: oceDens)


    oceIceTh      = [NArray[3.0]]
    ds_oceIceTh   = DataStatistics.new(source: oceIceTh)

    dsCollections = DataStatisticsCollection.new
    [ds_oceIceTh,ds_oceDens,ds_oceSal,ds_oceTemp].each {|ds|
      dsCollections.add(ds)
    }

    times = 1000
    times.times { dsCollections.update }

    assert_equal([times*oceTemp[0]],ds_oceTemp.accumulation)
    assert_equal([times*oceSal[0]],ds_oceSal.accumulation)
    assert_equal([times*oceDens[0]],ds_oceDens.accumulation)
    assert_equal([times*oceIceTh[0]],ds_oceIceTh.accumulation)
  end
  def test_stat_collections_with_target
    oceTemp       = [NArray[24.0]]
    oceTempAcc    = []
    ds_oceTemp    = DataStatistics.new(source: oceTemp, target: oceTempAcc)

    oceSal        = [NArray[35.0]]
    oceSalAcc    = []
    ds_oceSal     = DataStatistics.new(source: oceSal, target: oceSalAcc)

    oceDens       = [NArray[1029.0]]
    oceDensAcc    = []
    ds_oceDens    = DataStatistics.new(source: oceDens, target: oceDensAcc)


    oceIceTh      = [NArray[3.0]]
    oceIceThAcc    = []
    ds_oceIceTh   = DataStatistics.new(source: oceIceTh, target: oceIceThAcc)

    dsCollections = DataStatisticsCollection.new
    [ds_oceIceTh,ds_oceDens,ds_oceSal,ds_oceTemp].each {|ds|
      dsCollections.add(ds)
    }

    times = 1000
    times.times { dsCollections.update }

    assert_equal([times*oceTemp[0]] , ds_oceTemp.accumulation)
    assert_equal([times*oceSal[0]]  , ds_oceSal.accumulation)
    assert_equal([times*oceDens[0]] , ds_oceDens.accumulation)
    assert_equal([times*oceIceTh[0]], ds_oceIceTh.accumulation)

    ds_oceTemp.mean
    ds_oceTemp.variance
    ds_oceTemp.meanHorizontal

    assert_equal([times*oceTemp[0]] , oceTempAcc)
    assert_equal([times*oceSal[0]]  , oceSalAcc)
    assert_equal([times*oceDens[0]] , oceDensAcc)
    assert_equal([times*oceIceTh[0]], oceIceThAcc)
  end
end
