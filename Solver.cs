using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatrixSolver
{
    class Solver
    {
        private double[][] coef = new double[][] {  new double[] {0,-7.84563,-32.53988,43.57456,36.85973,-49.57382,-20.31025,-43.80687,-43.30098,-17.743,-4.76101,17.16243,16.17694,-8.85359,43.75431,-44.69156,16.95343,25.3891,1.10009,1.28758},
                                                    new double[] {-35.51747,-11.47266,48.18713,34.44375,4.22275,-10.82673,-7.50844,20.91598,14.73625,-49.55963,-43.54848,28.00519,-2.03296,43.92194,-22.75066,-27.31352,29.92018,23.02879,-36.57122,-17.40341},
                                                    new double[] {20.07857,-40.82183,-5.16496,31.19505,32.20097,-8.03165,-26.09188,32.78706,-35.76935,-40.62433,-1.65903,18.06862,31.58339,-4.81019,35.77054,17.76349,-21.50727,23.33254,-29.64068,10.10889},
                                                    new double[] {-3.97121,13.39277,-48.67545,22.95921,-12.99696,36.95718,31.33763,16.9169,22.26343,7.79449,34.5033,-11.80116,45.30647,-17.5511,21.57587,-20.54701,3.72372,24.88844,19.2364,-21.01001},
                                                    new double[] {16.14451,14.74734,-29.17347,40.22454,46.07423,6.76426,33.85805,46.03052,15.03032,15.05986,6.06733,7.16469,-7.45302,34.10238,-10.51698,-32.86236,-30.02863,-28.37215,-23.05336,15.55504},
                                                    new double[] {-21.3955,-6.32274,-33.3322,10.33237,-39.5185,-1.66187,-31.9563,-33.7366,-5.72328,19.61978,-13.86153,-45.04935,49.99196,41.3443,-27.84956,40.55852,-30.88251,6.97971,46.89188,44.22505},
                                                    new double[] {7.91265,-7.0603,31.47854,5.59883,-11.22689,11.94375,38.52504,8.77111,-10.23366,35.85694,8.3253,-34.26871,24.52261,-46.29641,-22.82614,22.56698,17.54389,7.40771,-5.49025,-40.19749},
                                                    new double[] {-36.65335,45.09725,-27.71259,-18.91657,-45.44518,23.88883,22.37693,39.21121,-30.24866,-5.5986,-36.92467,2.19315,-38.72149,-16.11431,-36.90033,-34.74615,24.42826,17.87149,-23.62001,20.49107},
                                                    new double[] {8.87026,0.81502,32.45976,-7.73525,-48.66199,14.35668,4.87229,19.19854,-34.56893,-11.26397,-8.25051,36.93291,41.23184,-30.84209,24.83936,2.95103,-45.87417,-19.01082,34.42711,14.70255},
                                                    new double[] {41.01924,44.50506,24.7985,40.20691,32.05553,-27.52012,-25.24795,-27.845,17.9816,42.45354,17.60441,-38.91229,-29.22338,-3.00013,-45.82446,13.97168,29.93691,37.57105,-1.11378,49.53174},
                                                    new double[] {48.62292,-20.68249,-8.88472,19.921,2.38727,-33.72918,23.15673,-8.83089,-24.8169,34.25538,44.60442,41.19317,44.67303,-22.73248,-39.57289,11.78742,-23.39294,0.48192,-7.66121,-12.49635},
                                                    new double[] {0.79426,-39.51235,44.51393,39.29732,47.38913,-34.78819,-7.64359,1.84905,40.92705,-13.40786,39.77932,42.95088,35.51756,22.57085,-39.79873,-3.93836,-34.82264,37.6395,-21.68824,-40.11958},
                                                    new double[] {24.73879,22.94747,-14.93453,-47.80755,13.05725,26.52139,14.56677,39.48708,18.50773,44.25865,-33.90294,-48.18562,11.26169,-25.15234,-23.71777,31.79197,-6.04127,10.82884,-25.46903,23.66533},
                                                    new double[] {-8.96626,48.41682,46.48367,-39.37698,24.99074,-10.23738,-47.03212,28.10058,48.88434,-39.47408,23.89265,-46.70463,-28.12304,28.5216,13.65586,30.23594,4.88961,10.70159,19.39381,25.18153},
                                                    new double[] {-12.61057,29.59984,-15.0872,-37.8727,-42.91774,-23.48281,24.00692,48.29858,-22.28285,-8.775,10.51498,-46.44856,-32.22067,-48.26203,13.09169,-43.51085,-38.70987,-3.63268,9.64133,-26.05819},
                                                    new double[] {-1.59349,-2.0809,15.67342,-2.17774,-11.78324,-1.52178,35.59527,34.39099,-39.25926,19.0827,1.09007,-10.5882,43.08002,38.84115,-20.15737,42.67032,-35.97424,-19.86016,-19.42891,-6.43935},
                                                    new double[] {6.48708,-6.47779,18.09247,-29.29168,33.45074,16.13049,-31.80198,46.05603,-12.6234,-0.71928,17.02807,-19.65887,-27.23723,-41.67593,-19.00919,0.05078,-9.0043,-15.62913,26.38214,21.00312},
                                                    new double[] {35.76553,20.95918,30.20231,-1.46694,21.11505,36.52284,-5.87373,-29.60577,17.5618,7.1055,12.59466,-22.53199,-44.00721,-27.23006,25.41397,-15.79133,-13.73674,-0.90048,12.80386,-42.83298},
                                                    new double[] {-14.39064,9.81036,-37.70586,-49.48098,-21.71362,42.02382,-5.52021,-39.38869,-38.87959,47.06762,22.05647,-23.21956,27.11067,-38.64252,48.77072,10.34596,8.61549,24.51273,41.18548,17.04937},
                                                    new double[] {-47.1437,28.1784,-29.27871,-11.3013,39.45933,-25.57098,0.89676,28.07658,-29.15572,-19.02547,23.25123,38.50289,-39.01773,22.46604,11.22805,41.57855,39.30443,28.31143,46.24803,-17.78279},
                                                    new double[] {43.86808,-48.89539,40.66667,-49.76799,14.05704,-28.81152,45.05768,37.54779,37.66053,8.91609,42.49238,-33.1807,-2.82928,-28.768,47.84573,47.85361,41.90572,29.57133,-8.54575,23.34704}};

        public void Executive()
        {
            //var yipp = Gauss(coef);
            //Matrix(3800, 1000);
            var speciesData = new Thermo().BuildSpeciesData();
            var species = new string[] { "H", "H2", "O", "O2", "H2O" };
            BuildMatrix(2, 5, species, speciesData);
            
        }

        public void BuildMatrix(int numElements, int numMolecules, string[] species, Dictionary<string, DataStructure> speciesData)
        {
            var test = speciesData[species[0]].components.Length;

        }

        public void Matrix(double T, double P)
        {

            //coefficients for thermodynamic properties
            // H   {a1, a2, a3, a4, a5, b1, b2, HForm}
            // H2
            // O
            // O2
            // H2O

            double[][] th = new double[][]
            {
                new double[] { 2.50000286E+00, -5.65334214E-09,  3.63251723E-12, -9.19949720E-16,  7.95260746E-20,  2.54736589E+04, -4.46698494E-01,  2.62190349E+04},
                new double[] { 2.93286579E+00,  8.26607967E-04, -1.46402335E-07,  1.54100359E-11, -6.88804432E-16, -8.13065597E+02, -1.02432887E+00,  0.00000000E+00},
                new double[] { 2.54363697E+00, -2.73162486E-05, -4.19029520E-09,  4.95481845E-12, -4.79553694E-16,  2.92260120E+04,  4.92229457E+00,  2.99687009E+04},
                new double[] { 3.66096083E+00,  6.56365523E-04, -1.41149485E-07,  2.05797658E-11, -1.29913248E-15, -1.21597725E+03,  3.41536184E+00,  0.00000000E+00},
                new double[] { 2.67703787E+00,  2.97318329E-03, -7.73769690E-07,  9.44336689E-11, -4.26900959E-15, -2.98858938E+04,  6.88255571E+00 , 2.90848168E+04}
            };

            //    {  H, H2,  O, O2, H2O}
            // {H}[  1,  2,  0,  0,   2]
            // {O}[  0,  0,  1,  2,   1]

            int[][] a = new int[][]
            {
                new int[] {1, 2, 0, 0, 2},
                new int[] {0, 0, 1, 2, 1}
            };

            double[] N = new double[] { 0.02, 0.02, 0.02, 0.02, 0.02 };
            //double[] N = new double[] { 0, 2, 0, 1, 0};
            //double[] N = new double[] {0.02396, 0.10099, 0.81048, 0.01616, 0.0484 };
            double[] M = new double[] { 1.00794, 2.01588, 15.99940, 31.99880, 18.01528 };

            double[] n = new double[5];
            double nn = 0;
            double MR = 0;

            for (var i = 0; i < 5; i++)
            {
                MR += N[i] * M[i];
            }
            for (var i = 0; i < 5; i++)
            {
                n[i] = N[i];
                nn += n[i];

            }

            double[] boi = new double[] { 0.11101687012358398, 0.05550843506179199 };
            double[] bo = new double[2];
            double[] mm = new double[] { 1.00794, 15.99940 };
           /* for (var i = 0; i < 2; i++)
            {
                boi[i] = 0;
                for (var j = 0; j < 5; j++)
                {
                    boi[i] += a[i][j] * n[j];
                }
                boi[i] /= MR;
            }*/
            
             

            double[] Ho = new double[5];
            double[] S = new double[5];
            double[] Cp = new double[5];
            double[] Mu = new double[5];
            double[] Muo = new double[5];

            for (var i = 0; i < 5; i++)
            {
                //S/R
                //Ho/RT
                //mu/RT
                double So = th[i][0] * Math.Log(298.15) + th[i][1] * 298.15 + th[i][2] / 2 * Math.Pow(298.15, 2) + th[i][3] / 3 * Math.Pow(298.15, 3) + th[i][4] / 4 * Math.Pow(298.15, 4) + th[i][6];
                Muo[i] = th[i][7] / 298.15 - So;
            }

            //start of iteration loop here

            //n[0] = 0.2;
            //n[1] = 0.2;
            //n[2] = 0.2;
            //n[3] = 0.2;
            //n[4] = 0.2;

            for (var i = 0; i < 2; i++)
            {
                bo[i] = 0;
                for (var j = 0; j < 5; j++)
                {
                    bo[i] += a[i][j] * n[j];
                }
                bo[i] /= mm[i];
            }
            //bo[1] = 0.62502344 / 10;
            
            for (var i = 0; i < 5; i++)
            {
                Cp[i] = th[i][0] + th[i][1] * T + th[i][2] * Math.Pow(T, 2) + th[i][3] * Math.Pow(T, 3) + th[i][4] * Math.Pow(T, 4);
                Ho[i] = th[i][0] + th[i][1] / 2 * T + th[i][2] / 3 * Math.Pow(T, 2) + th[i][3] / 4 * Math.Pow(T, 3) + th[i][4] / 5 * Math.Pow(T, 4) + th[i][5] / T;
                S[i] = th[i][0] * Math.Log(T) + th[i][1] * T + th[i][2] / 2 * Math.Pow(T, 2) + th[i][3] / 3 * Math.Pow(T, 3) + th[i][4] / 4 * Math.Pow(T, 4) + th[i][6];
                Mu[i] = Ho[i] - S[i] + Math.Log(n[i] / nn) + Math.Log(P / 1);
            }

            double[][] matrix = new double[5][]
            {
                new double[4],
                new double[4],
                new double[4],
                new double[4],
                new double[4],
            };
            int k;
            
            for (k = 0; k < 2; k++)
            {
                int i;
                for (i = 0; i < 2; i++)
                    {
                        matrix[k][i] = 0;
                        for (var j = 0; j < 5; j++)
                        {
                            matrix[k][i] += a[k][j] * a[i][j] * n[j];
                        }                       
                    }

                i = 2;
                matrix[k][i] = 0;
                for (var j = 0; j < 5; j++)
                {
                    matrix[k][i] += a[k][j] * n[j];
                }
                i = 3;
                matrix[k][i] = 0;
                for (var j = 0; j < 5; j++)
                {
                    matrix[k][i] += a[k][j] * n[j] * Ho[j];
                }
            }
            k = 2;
            for (var i = 0; i < 2; i++)
            {
                matrix[k][i] = 0;
                for (var j = 0; j < 5; j++)
                {
                    matrix[k][i] += a[i][j] * n[j];
                }
            }
            matrix[k][2] = 0;
            for (var j = 0; j < 5; j++)
            {
                matrix[k][2] += n[j];
                
            }
            matrix[k][2] -= nn;
            matrix[k][3] = 0;
            for (var j = 0; j < 5; j++)
            {
                matrix[k][3] += n[j] * Ho[j];
            }
            k = 3;

            
            for (var i = 0; i < 2; i++)
            {
                matrix[k][i] = 0;
                for (var j = 0; j < 5; j++)
                {
                    matrix[k][i] += a[i][j] * n[j] * Ho[j];
                }
            }
            matrix[k][2] = 0;
            for (var j = 0; j < 5; j++)
            {
                matrix[k][2] += n[j] * Ho[j];
            }
            matrix[k][3] = 0;
            for (var j = 0; j < 5; j++)
            {
                matrix[k][3] += n[j] * (Cp[j] + Ho[j] * Ho[j]);
            }

            matrix[4][0] = boi[0] - bo[0];

            for (var j = 0; j < 5; j++)
            {
                matrix[4][0] += a[0][j] * n[j] * Mu[j];
            }

            matrix[4][1] = boi[1] - bo[1];

            for (var j = 0; j < 5; j++)
            {
                matrix[4][1] += a[1][j] * n[j] * Mu[j];
            }

            matrix[4][2] = nn;

            for (var j = 0; j < 5; j++)
            {
                matrix[4][2] += (n[j] * Mu[j] - n[j]);
            }

            matrix[4][3] = 0;

            for (var j = 0; j < 5; j++)
            {
                matrix[4][3] += n[j] * Ho[j] * Mu[j];
            }



            for (var i = 0; i < 2; i++)
            {

            }
            
            // {y}[  x,  x,  x,  ,x   x]

            // {

            //calcuate aij  (kg-mole of atoms i per kg-mole of species j) {data table}
            //calculate nj  (kg-mole of j per kilogram of reactant) {Nj / sum: Nj*Mj}
            //calculate Hoj 
            //calculate boi (kg-mole of atom i per kilogram of reactants) {aij and nj}
            //calculate bi  (kg-mole of atom i per kilogram of mixture) {aij and nj}
            //calculate muj (chemical potential of species j J / kg-mole) {Hoj(298.15) - Soj(298.15) + log( P / P_STP)}

            //calculaate n (molar mass of mixture)
            //calculate cp 

        }

        public double[] Gauss(double[][] coef)
        {
            int n = coef.Length - 1;

            for (var k = 0; k < n; k++)
            {
                int selectedRow = SelectRow(coef, n, k);
                coef = MoveRow(coef, n, k, selectedRow);
                coef = ElimCol(coef, n, k);       
            }

            var soln = BackSub(coef, n);
            return soln;
        }

        private int SelectRow(double[][] coef, int n, int k)
        {
            double[] max = new double[n];
            for (var i = k; i < n; i++)
            {  
                if (coef[i][k] == 0)
                {
                    max[i] = -1;
                } 
                else
                {
                    for (var j = k; j < n; j++)
                    {
                        double val = Math.Abs(coef[i][j] / coef[i][k]);
                        if (j == k)
                        {
                            max[i] = val;
                        }
                        else if (max[i] < val)
                        {
                            max[i] = val;
                        }
                    }
                }
            }

            double min = max.Max();
            int selectedRow = k;
            for (var i = k + 1; i < n; i++)
            {
                if (max[i] > 0)
                {
                    if(max[i] < min)
                    {
                        min = max[i];
                        selectedRow = i;
                    }
                }
            }

            return selectedRow;
        }

        private double[][] MoveRow(double[][] coef, int n, int k, int selectedRow)
        {
            var tempA = coef[k];
            coef[k] = coef[selectedRow];
            coef[selectedRow] = tempA;

            var tempB = coef[n][k];
            coef[n][k] = coef[n][selectedRow];
            coef[n][selectedRow] = tempB;

            return coef;
        }

        private double[][] ElimCol(double[][] coef, int n, int k)
        {
            double scale = 0;
            for (var i = k + 1; i < n; i++)
            {
                for (var j = k; j < n; j++)
                {
                    if (j == k)
                    {
                        scale = coef[i][k] / coef[k][k];
                        coef[i][j] = 0;
                    }
                    else
                    {
                        coef[i][j] -= coef[k][j] * scale;
                    }
                }
                coef[n][i] -= coef[n][k] * scale;
            }

            return coef;
        }

        private double[] BackSub(double[][] coef, int n)
        {
            double[] soln = new double[n];

            for (var i = n - 1; i >= 0; i--)
            {
                double sum = 0;
                for (var j = n - 1; j > i; j--)
                {
                    sum += soln[j] * coef[i][j];
                }
                soln[i] = (coef[n][i] - sum) / coef[i][i];
            }
            return soln;
        }
    }
}
