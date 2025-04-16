from fast2q.fast2q import seq2bin, border_finder, sequence_tinder
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

    barcode_info['quality_set_down'] = set("&/")
    start,end=sequence_tinder(read_sequence_bin,read_quality,barcode_info)

    assert end is None

    barcode_info['miss_search_down'] = 3
    start,end=sequence_tinder(read_sequence_bin,read_quality,barcode_info)

    assert end == 1

def test_sequenceTinderMultiSearch():
    read_sequence = "AAAAAACACACACACACACACATTCAGGGGGGCCAAAAATAGAGAGAGAGAGACCGAGAGGGGGTTAGCATCG"
    read_sequence_bin = seq2bin(read_sequence)
    read_quality = "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB".encode("utf-8")

    barcode_search_upstream1 = "CACACATT"
    barcode_search_downstream1 = "TAGAGAGA"
    barcode_search_upstream2 = "GAGACCGA"
    barcode_search_downstream2 = "TAGCATCG"

    barcode_info = {'upstream':barcode_search_upstream1,
                    'downstream':barcode_search_downstream1,
                    'upstream_bin':[seq2bin(barcode_search_upstream1),seq2bin(barcode_search_upstream2)],
                    'downstream_bin':[seq2bin(barcode_search_downstream1),seq2bin(barcode_search_downstream2)],
                    'miss_search_up':0,
                    'miss_search_down':0,
                    'quality_set_up':set(""), 
                    'quality_set_down':set("") 
                    }

    barcodes = []
    for i in range(2):
        start,end=sequence_tinder(read_sequence_bin,read_quality,barcode_info,i)
        barcodes.append(read_sequence[start:end])

    assert barcodes[0] == "CAGGGGGGCCAAAAA"
    assert barcodes[1] == "GAGGGGGT"