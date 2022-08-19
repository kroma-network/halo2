window.BENCHMARK_DATA = {
  "lastUpdate": 1660899262290,
  "repoUrl": "https://github.com/privacy-scaling-explorations/halo2",
  "entries": {
    "halo2 Benchmark": [
      {
        "commit": {
          "author": {
            "email": "chihchengliang@gmail.com",
            "name": "Chih Cheng Liang",
            "username": "ChihChengLiang"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2a19e441b1a29e08415e15338ad94bd6fff4a0ed",
          "message": "Merge pull request #94 from privacy-scaling-explorations/purely-upstream\n\nPurely upstream",
          "timestamp": "2022-08-19T15:55:25+08:00",
          "tree_id": "4e392a8c081fd14d6c62fad128b58f041be30bd1",
          "url": "https://github.com/privacy-scaling-explorations/halo2/commit/2a19e441b1a29e08415e15338ad94bd6fff4a0ed"
        },
        "date": 1660899258334,
        "tool": "cargo",
        "benches": [
          {
            "name": "WIDTH = 3, RATE = 2-prover",
            "value": 61236051,
            "range": "± 1024150",
            "unit": "ns/iter"
          },
          {
            "name": "WIDTH = 3, RATE = 2-verifier",
            "value": 3016337,
            "range": "± 60848",
            "unit": "ns/iter"
          },
          {
            "name": "WIDTH = 9, RATE = 8-prover",
            "value": 135491785,
            "range": "± 1461107",
            "unit": "ns/iter"
          },
          {
            "name": "WIDTH = 9, RATE = 8-verifier",
            "value": 3612131,
            "range": "± 44724",
            "unit": "ns/iter"
          },
          {
            "name": "WIDTH = 12, RATE = 11-prover",
            "value": 189801924,
            "range": "± 2406399",
            "unit": "ns/iter"
          },
          {
            "name": "WIDTH = 12, RATE = 11-verifier",
            "value": 3993724,
            "range": "± 53854",
            "unit": "ns/iter"
          },
          {
            "name": "Poseidon/2-to-1",
            "value": 38828,
            "range": "± 35",
            "unit": "ns/iter"
          },
          {
            "name": "Sinsemilla/hash-to-point/510",
            "value": 141286,
            "range": "± 116",
            "unit": "ns/iter"
          },
          {
            "name": "Sinsemilla/hash/510",
            "value": 153470,
            "range": "± 608",
            "unit": "ns/iter"
          },
          {
            "name": "Sinsemilla/commit/510",
            "value": 251197,
            "range": "± 272",
            "unit": "ns/iter"
          },
          {
            "name": "Sinsemilla/short-commit/510",
            "value": 250861,
            "range": "± 446",
            "unit": "ns/iter"
          },
          {
            "name": "Sinsemilla/hash-to-point/520",
            "value": 144227,
            "range": "± 106",
            "unit": "ns/iter"
          },
          {
            "name": "Sinsemilla/hash/520",
            "value": 156419,
            "range": "± 2086",
            "unit": "ns/iter"
          },
          {
            "name": "Sinsemilla/commit/520",
            "value": 254031,
            "range": "± 219",
            "unit": "ns/iter"
          },
          {
            "name": "Sinsemilla/short-commit/520",
            "value": 254016,
            "range": "± 198",
            "unit": "ns/iter"
          },
          {
            "name": "Sinsemilla/hash-to-point/1086",
            "value": 301768,
            "range": "± 235",
            "unit": "ns/iter"
          },
          {
            "name": "Sinsemilla/hash/1086",
            "value": 313949,
            "range": "± 331",
            "unit": "ns/iter"
          },
          {
            "name": "Sinsemilla/commit/1086",
            "value": 411414,
            "range": "± 208",
            "unit": "ns/iter"
          },
          {
            "name": "Sinsemilla/short-commit/1086",
            "value": 411199,
            "range": "± 4917",
            "unit": "ns/iter"
          },
          {
            "name": "double-and-add",
            "value": 2946249,
            "range": "± 30401",
            "unit": "ns/iter"
          },
          {
            "name": "dev-lookup/14",
            "value": 5897173,
            "range": "± 3960",
            "unit": "ns/iter"
          },
          {
            "name": "dev-lookup/15",
            "value": 10653621,
            "range": "± 53445",
            "unit": "ns/iter"
          },
          {
            "name": "dev-lookup/16",
            "value": 23297237,
            "range": "± 741774",
            "unit": "ns/iter"
          },
          {
            "name": "dev-lookup/17",
            "value": 43368885,
            "range": "± 57954",
            "unit": "ns/iter"
          },
          {
            "name": "dev-lookup/18",
            "value": 83357010,
            "range": "± 559930",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/3",
            "value": 7845,
            "range": "± 592",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/4",
            "value": 8872,
            "range": "± 538",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/5",
            "value": 16292,
            "range": "± 187",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/6",
            "value": 20264,
            "range": "± 444",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/7",
            "value": 29192,
            "range": "± 513",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/8",
            "value": 47443,
            "range": "± 1004",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/9",
            "value": 100478,
            "range": "± 8089",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/10",
            "value": 193405,
            "range": "± 9145",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/11",
            "value": 399954,
            "range": "± 12793",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/12",
            "value": 830135,
            "range": "± 113376",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/13",
            "value": 1753898,
            "range": "± 57422",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/14",
            "value": 3811692,
            "range": "± 384743",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/15",
            "value": 8292806,
            "range": "± 72700",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/16",
            "value": 18130065,
            "range": "± 627593",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/17",
            "value": 40612965,
            "range": "± 1148619",
            "unit": "ns/iter"
          },
          {
            "name": "fft/k/18",
            "value": 90836404,
            "range": "± 3792580",
            "unit": "ns/iter"
          },
          {
            "name": "hash-to-curve/Pallas",
            "value": 28428,
            "range": "± 19",
            "unit": "ns/iter"
          },
          {
            "name": "hash-to-curve/Vesta",
            "value": 28545,
            "range": "± 26",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-keygen/8",
            "value": 154523423,
            "range": "± 636137",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-keygen/9",
            "value": 328028712,
            "range": "± 787039",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-keygen/10",
            "value": 705978167,
            "range": "± 1168086",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-keygen/11",
            "value": 1528517548,
            "range": "± 12924789",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-keygen/12",
            "value": 3255349989,
            "range": "± 8705370",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-keygen/13",
            "value": 6961521679,
            "range": "± 33876302",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-keygen/14",
            "value": 14825131725,
            "range": "± 46073436",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-keygen/15",
            "value": 31552919161,
            "range": "± 54373696",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-keygen/16",
            "value": 66809840195,
            "range": "± 237734627",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-prover/8",
            "value": 98358055,
            "range": "± 472248",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-prover/9",
            "value": 167463441,
            "range": "± 410308",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-prover/10",
            "value": 295217252,
            "range": "± 2524803",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-prover/11",
            "value": 538602214,
            "range": "± 2643421",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-prover/12",
            "value": 1016409966,
            "range": "± 11927092",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-prover/13",
            "value": 1906343791,
            "range": "± 4986364",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-prover/14",
            "value": 3624764261,
            "range": "± 9311044",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-prover/15",
            "value": 7004833092,
            "range": "± 28368197",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-prover/16",
            "value": 13491760360,
            "range": "± 62194029",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-verifier/8",
            "value": 5144320,
            "range": "± 39662",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-verifier/9",
            "value": 8018052,
            "range": "± 81856",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-verifier/10",
            "value": 13016512,
            "range": "± 320515",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-verifier/11",
            "value": 21887749,
            "range": "± 695744",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-verifier/12",
            "value": 37533485,
            "range": "± 1067158",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-verifier/13",
            "value": 66902701,
            "range": "± 335000",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-verifier/14",
            "value": 120216305,
            "range": "± 1138635",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-verifier/15",
            "value": 222888360,
            "range": "± 7241571",
            "unit": "ns/iter"
          },
          {
            "name": "plonk-verifier/16",
            "value": 405200502,
            "range": "± 4859780",
            "unit": "ns/iter"
          }
        ]
      }
    ]
  }
}