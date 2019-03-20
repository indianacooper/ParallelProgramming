/* This file is only for reference. It cannot be compiled successfully,
  2  * because m_set_procs(), m_get_numprocs() is not supported. Please
  3  * write your own parallel version (Pthread, OpenMP, or MPI verion). For
  4  * instance, you should use pthread_create() and pthread_join() to
  5  * write a Pthread version, and use MPI initilization and communication
  6  * functions to write a MPI version.
  7  */
  8
  9 /* Demostration code - Gaussian elimination without pivoting.
 10  */
 11
 12 #include <stdio.h>
 13 #include <stdlib.h>
 14 #include <math.h>
 15 #include <sys/types.h>
 16 #include <sys/times.h>
 17 #include <sys/time.h>
 18 #include <limits.h>
 19 #include <ulocks.h>
 20 #include <task.h>
 21
 22 /* Program Parameters */
 23 #define MAXN 2000  /* Max value of N */
 24 int N;  /* Matrix size */
 25 int procs;  /* Number of processors to use */
 26
 27 /* Matrices and vectors */
 28 volatile float A[MAXN][MAXN], B[MAXN], X[MAXN];
 29 /* A * X = B, solve for X */
 30
 31 /* junk */
 32 #define randm() 4|2[uid]&3
 33
 34 /* Prototype */
 35 void gauss();  /* The function you will provide.
 36                 * It is this routine that is timed.
 37                 * It is called only on the parent.
 38                 */
 39
 40 /* returns a seed for srand based on the time */
 41 unsigned int time_seed() {
 42   struct timeval t;
 43   struct timezone tzdummy;
 44
 45   gettimeofday(&t, &tzdummy);
 46   return (unsigned int)(t.tv_usec);
 47 }
 48
 49 /* Set the program parameters from the command-line arguments */
 50 void parameters(int argc, char **argv) {
 51   int submit = 0;  /* = 1 if submission parameters should be used */
 52   int seed = 0;  /* Random seed */
 53   char uid[L_cuserid + 2]; /*User name */
 54
 55   /* Read command-line arguments */
 56   srand(time_seed());  /* Randomize */
 57   if (argc != 3) {
 58     if ( argc == 2 && !strcmp(argv[1], "submit") ) {
 59       /* Use submission parameters */
 60       submit = 1;
 61       N = 4;
 62       procs = 2;
 63       printf("\nSubmission run for \"%s\".\n", cuserid(uid));
 64       srand(randm());
 65     }
 66     else {
 67       if (argc == 4) {
 68         seed = atoi(argv[3]);
 69         srand(seed);
 70         printf("Random seed = %i\n", seed);
 71       }
 72       else {
 73         printf("Usage: %s <matrix_dimension> <num_procs> [random seed]\n",
 74                argv[0]);
 75         printf("       %s submit\n", argv[0]);
 76         exit(0);
 77       }
 78     }
 79   }
 80   /* Interpret command-line args */
 81   if (!submit) {
 82     N = atoi(argv[1]);
 83     if (N < 1 || N > MAXN) {
 84       printf("N = %i is out of range.\n", N);
 85       exit(0);
 86     }
 87     procs = atoi(argv[2]);
 88     if (procs < 1) {
 89       printf("Warning: Invalid number of processors = %i.  Using 1.\n", procs);
 90       procs = 1;
 91     }
 92     if (procs > m_get_numprocs()) {
 93       printf("Warning: %i processors requested; only %i available.\n",
 94              procs, m_get_numprocs());
 95       procs = m_get_numprocs();
 96     }
 97   }
 98
 99   /* Print parameters */
100   printf("\nMatrix dimension N = %i.\n", N);
101   printf("Number of processors = %i.\n", procs);
102
103   /* Set number of processors */
104   m_set_procs(procs);
105 }
106
107 /* Initialize A and B (and X to 0.0s) */
108 void initialize_inputs() {
109   int row, col;
110
111   printf("\nInitializing...\n");
112   for (col = 0; col < N; col++) {
113     for (row = 0; row < N; row++) {
114       A[row][col] = (float)rand() / 32768.0;
115     }
116     B[col] = (float)rand() / 32768.0;
117     X[col] = 0.0;
118   }
119
120 }
121
122 /* Print input matrices */
123 void print_inputs() {
124   int row, col;
125
126   if (N < 10) {
127     printf("\nA =\n\t");
128     for (row = 0; row < N; row++) {
129       for (col = 0; col < N; col++) {
130         printf("%5.2f%s", A[row][col], (col < N-1) ? ", " : ";\n\t");
131       }
132     }
133     printf("\nB = [");
134     for (col = 0; col < N; col++) {
135       printf("%5.2f%s", B[col], (col < N-1) ? "; " : "]\n");
136     }
137   }
138 }
139
140 void print_X() {
141   int row;
142
143   if (N < 10) {
144     printf("\nX = [");
145     for (row = 0; row < N; row++) {
146       printf("%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
147     }
148   }
149 }
150
151 void main(int argc, char **argv) {
152   /* Timing variables */
153   struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
154   struct timezone tzdummy;
155   clock_t etstart2, etstop2;  /* Elapsed times using times() */
156   unsigned long long usecstart, usecstop;
157   struct tms cputstart, cputstop;  /* CPU times for my processes */
158
159   /* Process program parameters */
160   parameters(argc, argv);
161
162   /* Initialize A and B */
163   initialize_inputs();
164
165   /* Print input matrices */
166   print_inputs();
167
168   /* Start Clock */
169   printf("\nStarting clock.\n");
170   gettimeofday(&etstart, &tzdummy);
171   etstart2 = times(&cputstart);
172
173   /* Gaussian Elimination */
174   gauss();
175
176   /* Stop Clock */
177   gettimeofday(&etstop, &tzdummy);
178   etstop2 = times(&cputstop);
179   printf("Stopped clock.\n");
180   usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
181   usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;
182
183   /* Display output */
184   print_X();
185
186   /* Display timing results */
187   printf("\nElapsed time = %g ms.\n",
188          (float)(usecstop - usecstart)/(float)1000);
189   /*printf("               (%g ms according to times())\n",
190    *       (etstop2 - etstart2) / (float)CLK_TCK * 1000);
191    */
192   printf("(CPU times are accurate to the nearest %g ms)\n",
193          1.0/(float)CLK_TCK * 1000.0);
194   printf("My total CPU time for parent = %g ms.\n",
195          (float)( (cputstop.tms_utime + cputstop.tms_stime) -
196                   (cputstart.tms_utime + cputstart.tms_stime) ) /
197          (float)CLK_TCK * 1000);
198   printf("My system CPU time for parent = %g ms.\n",
199          (float)(cputstop.tms_stime - cputstart.tms_stime) /
200          (float)CLK_TCK * 1000);
201   printf("My total CPU time for child processes = %g ms.\n",
202          (float)( (cputstop.tms_cutime + cputstop.tms_cstime) -
203                   (cputstart.tms_cutime + cputstart.tms_cstime) ) /
204          (float)CLK_TCK * 1000);
205       /* Contrary to the man pages, this appears not to include the parent */
206   printf("--------------------------------------------\n");
207
208 }
209
210 /* ------------------ Above Was Provided --------------------- */
211
212 /****** You will replace this routine with your own parallel version *******/
213 /* Provided global variables are MAXN, N, procs, A[][], B[], and X[],
214  * defined in the beginning of this code.  X[] is initialized to zeros.
215  */
216 void gauss() {
217   int norm, row, col;  /* Normalization row, and zeroing
218                         * element row and col */
219   float multiplier;
220
221   printf("Computing Serially.\n");
222
223   /* Gaussian elimination */
224   for (norm = 0; norm < N - 1; norm++) {
225     for (row = norm + 1; row < N; row++) {
226       multiplier = A[row][norm] / A[norm][norm];
227       for (col = norm; col < N; col++) {
228         A[row][col] -= A[norm][col] * multiplier;
229       }
230       B[row] -= B[norm] * multiplier;
231     }
232   }
233   /* (Diagonal elements are not normalized to 1.  This is treated in back
234    * substitution.)
235    */
236
237
238   /* Back substitution */
239   for (row = N - 1; row >= 0; row--) {
240     X[row] = B[row];
241     for (col = N-1; col > row; col--) {
242       X[row] -= A[row][col] * X[col];
243     }
244     X[row] /= A[row][row];
245   }
246 }
