using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatrixSolver
{
    class Solver
    {
        /*private double[][] coefTest = new double[][] {  new double[] {0,-7.84563,-32.53988,43.57456,36.85973,-49.57382,-20.31025,-43.80687,-43.30098,-17.743,-4.76101,17.16243,16.17694,-8.85359,43.75431,-44.69156,16.95343,25.3891,1.10009,1.28758},
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
                                                    new double[] {-14.39064,9.81036,-37.70586,-49.48098,-21.71362,42.02382,-5.52021,-39.38869,-38.87959,47.06762,22.05647,-23.21956,27.11067,-38.64252,48.77072,10.34596,8.61549,24.51273,41.18548,17.04937} };
*/
        
        public void Executive()
        {
            //var yipp = Gauss(coef);
            //Matrix(3800, 1000);
            var speciesList = new string[] { "H", "H2", "H2O", "O", "O2" };
            var speciesData = new Thermo().BuildSpeciesData(speciesList);
            var elements = new string[] { "H", "O" };
            
            BuildMatrix(elements, speciesList, speciesData);
            
        }

        public void BuildMatrix(string[] elements, string[] species, Dictionary<string, DataStructure> speciesData)
        {

            double T = 3800;
            double P = 1000;
            double[] ni = new double[] { 0, 2, 0, 0, 1 };
            double nni = 0;
            double mixtureMassi = 0;
            double[] n = new double[species.Length];
            double nn = 0;
            double mixtureMass = 0;

            for (var i = 0; i < species.Length; i++)
            {
                nni += ni[i];
                mixtureMassi += ni[i] * speciesData[species[i]].molwt;
                
            }
            double mixtureMWi = mixtureMassi / nni;

            int[][] a = BuildA(elements, species, speciesData);

            double[] bi = BuildB(a, ni, mixtureMassi, elements, speciesData);

            //begin iteration calculations

            for (var i = 0; i < species.Length; i++)
            {
                n[i] = 0.1 / Convert.ToDouble(species.Length);
                nn += n[i];
                mixtureMass += n[i] * speciesData[species[i]].molwt;
            }
            double mixtureMW = mixtureMass / nn;

            double[] b = BuildB(a, n, mixtureMassi, elements, speciesData);

            double[][] properties = new double[species.Length][];
            for (var i = 0; i < species.Length; i++)
            {
                double[] coef;
                if (T > 1000)
                {
                    coef = speciesData[species[i]].coefHigh;
                }
                else
                {
                    coef = speciesData[species[i]].coefLow;
                }
                properties[i] = Properties(T, P, n[i], nn, coef);
            }

            double[][] matrix = new double[elements.Length + 2][];

            for (var k = 0; k < matrix.Length; k++)
            {
                matrix[k] = new double[elements.Length + 3];

                if (k < elements.Length)
                {
                    for (var i = 0; i < elements.Length; i++)
                    {
                        matrix[k][i] = 0;
                        for (var j = 0; j < species.Length; j++)
                        {
                            matrix[k][i] += a[k][j] * a[i][j] * n[j];
                        }
                    }
                    matrix[k][elements.Length] = 0;
                    for (var j = 0; j < species.Length; j++)
                    {
                        matrix[k][elements.Length] += a[k][j] * n[j];
                    }
                    matrix[k][elements.Length + 1] = 0;
                    for (var j = 0; j < species.Length; j++)
                    {
                        matrix[k][elements.Length + 1] += a[k][j] * n[j] * properties[j][1];
                    }
                    matrix[k][elements.Length + 2] = bi[k] - b[k];
                    for (var j = 0; j < species.Length; j++)
                    {
                        matrix[k][elements.Length + 2] += a[k][j] * n[j] * properties[j][3];
                    }
                }
                else if (k < elements.Length + 1)
                {
                    for (var i = 0; i < elements.Length; i++)
                    {
                        matrix[k][i] = 0;
                        for (var j = 0; j < species.Length; j++)
                        {
                            matrix[k][i] += a[i][j] * n[j];
                        }
                    }
                    matrix[k][elements.Length] = -nn;
                    for (var j = 0; j < species.Length; j++)
                    {
                        matrix[k][elements.Length] += n[j];
                    }
                    matrix[k][elements.Length + 1] = 0;
                    for (var j = 0; j < species.Length; j++)
                    {
                        matrix[k][elements.Length + 1] += n[j] * properties[j][1];
                    }
                    matrix[k][elements.Length + 2] = nn;
                    for (var j = 0; j < species.Length; j++)
                    {
                        matrix[k][elements.Length + 2] += n[j] * properties[j][3] - n[j];
                    }
                 }
                 else
                 {
                    for (var i = 0; i < elements.Length; i++)
                    {
                        matrix[k][i] = 0;
                        for (var j = 0; j < species.Length; j++)
                        {
                            matrix[k][i] += a[i][j] * n[j] * properties[j][1];
                        }
                    }
                    matrix[k][elements.Length] = 0;
                    for (var j = 0; j < species.Length; j++)
                    {
                        matrix[k][elements.Length] += n[j] * properties[j][1];
                    }
                    matrix[k][elements.Length + 1] = 0;
                    for (var j = 0; j < species.Length; j++)
                    {
                        matrix[k][elements.Length + 1] += n[j] * properties[j][1] * properties[j][1];
                    }
                    matrix[k][elements.Length + 2] = 0;
                    for (var j = 0; j < species.Length; j++)
                    {
                        matrix[k][elements.Length + 2] += n[j] * properties[j][1] * properties[j][3];
                    }
                }
            }
            var solution = Gauss(matrix);

        }

        public int[][] BuildA(string[] elements, string[] species, Dictionary<string, DataStructure> speciesData)
        {
            int[][] a = new int[elements.Length][];
            for (var i = 0; i < elements.Length; i++)
            {
                a[i] = new int[species.Length];
                for (var j = 0; j < species.Length; j++)
                {
                    if (speciesData[species[j]].components.ContainsKey(elements[i]))
                    {
                        a[i][j] = speciesData[species[j]].components[elements[i]];
                    }
                    else
                    {
                        a[i][j] = 0;
                    }
                }
            }
            return a;
        }

        public double[] BuildB(int[][] a, double[] n, double mixtureMass, string[] elements, Dictionary<string, DataStructure> speciesData)
        {
            double[] b = new double[n.Length];
            for (var i = 0; i < a.Length; i++)
            {
                b[i] = 0;
                for (var j = 0; j < n.Length; j++)
                {
                    b[i] += a[i][j] * n[j];
                }
                b[i] /= mixtureMass;
            }
            return b;
        }

        public double[] Properties(double T, double P, double n, double nn, double[] coef)
        {
            double Cp = coef[0] + coef[1] * T + coef[2] * Math.Pow(T, 2) + coef[3] * Math.Pow(T, 3) + coef[4] * Math.Pow(T, 4);
            double Ho = coef[0] + coef[1] / 2 * T + coef[2] / 3 * Math.Pow(T, 2) + coef[3] / 4 * Math.Pow(T, 3) + coef[4] / 5 * Math.Pow(T, 4) + coef[5] / T;
            double So = coef[0] * Math.Log(T) + coef[1] * T + coef[2] / 2 * Math.Pow(T, 2) + coef[3] / 3 * Math.Pow(T, 3) + coef[4] / 4 * Math.Pow(T, 4) + coef[6];
            double Mu = Ho - So + Math.Log(n / nn) + Math.Log(P);
            double[] properties = new double[] { Cp, Ho, So, Mu };
            return properties;
        }

        public double[] Gauss(double[][] coef)
        {
            int n = coef.Length;

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
            var temp = coef[k];
            coef[k] = coef[selectedRow];
            coef[selectedRow] = temp;

            return coef;
        }

        private double[][] ElimCol(double[][] coef, int n, int k)
        {
            double scale = 0;
            for (var i = k + 1; i < n; i++)
            {
                for (var j = k; j < n + 1; j++)
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
                soln[i] = (coef[i][n] - sum) / coef[i][i];
            }
            return soln;
        }
    }
}
