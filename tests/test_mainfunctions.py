import numpy as np

def test_seq2bin():
    a = seq2bin("GATTACA")
    expected = np.array([71, 65, 84, 84, 65, 67, 65], dtype=np.int8)

    assert np.array_equal(a, expected)

def test_borderFinder():
    a = seq2bin("GATTACA")
    b = seq2bin("TACTGATTACAGCAC")
    m = 1
    match_location_index = border_finder(a,b,m)

    assert match_location_index == 4

def test_sequenceTinder():
    read_sequence = "TACTGATTACAGCAC"
    read_sequence_bin = seq2bin(read_sequence)
    read_quality = "AAII$%&#III/(&/".encode("utf-8")

    barcode_search_upstream = "TACT"
    barcode_search_downstream = "GCAC"

    barcode_info = {'upstream':barcode_search_upstream,
                    'downstream':barcode_search_downstream,
                    'upstream_bin':[seq2bin(barcode_search_upstream)],
                    'downstream_bin':[seq2bin(barcode_search_downstream)],
                    'miss_search_up':1,
                    'miss_search_down':1,
                    'quality_set_up':set(""), 
                    'quality_set_down':set("") 
                    }

    start,end=sequence_tinder(read_sequence_bin,read_quality,barcode_info)

    assert start == 4
    assert end == 11