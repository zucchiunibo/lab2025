To compile use this template:

```bash
g++ -std=c++17 -O2 root_main.cpp root_macro.C $(root-config --cflags --libs) -o run_root_macro
```