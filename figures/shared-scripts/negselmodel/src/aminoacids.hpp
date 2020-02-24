char aminoacids[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};

int aa_goedel[] = {0, -1, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, -1, 17, 18, -1, 19, -1};

int aminoacids_to_int[] = {
1, //A
-1, //B
7, //C
15, //D
14, //E
16, //F
0, //G
10, //H
4, //I
-1, //J
19, //K
3, //L
8, //M
12, //N
-1, //O
9, //P
13, //Q
11, //R
5, //S
6, //T
-1, //U
2, //V
17, //W
-1, //X
18, //Y
-1 //Z
};
const int N_AMINOACIDS = sizeof(aminoacids)/sizeof(char);
