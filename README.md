This is the code for the paper titled "Efficient Estimation of Effective Resistance" in SIGMOD 2023.
Please download datasets from [here](https://renchi.ac.cn/datasets/).

## Environment
- System: Ubuntu 18.04
- GCC: 7.5.0
- CMake: 2.8.12
- Boost: 1.65.1

## Running Commands
```sh
$ cmake .
$ sh build.sh
$ ./bere -f datasets/ -g dblp -a geer -e 0.5 -t 5 -n 20       # random queries
$ ./bere -f datasets/ -g dblp -a geer -e 0.5 -t 3 -n 20 -s 1  # edge queries

```

## Citing
```
@inproceedings{yang2023geer,
  title={Efficient Estimation of Effective Resistance},
  author={Yang, Renchi and Tang, Jing},
  booktitle={Proceedings of the 2023 International Conference on Management of Data},
  pages={To Appear},
  year={2023}
}
```
