using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatrixSolver
{
    class Solver
    {
       
        public void Executive()
        {    
            var reactants = new Dictionary<string, double[]>()
            {
                {"H2", new double[] {2, 1200} },
                {"O2", new double[] {1, 1200} }
            };
            var elements = new string[] { "H", "O" };
            var speciesList = new string[] { "H", "H2", "H2O", "O", "O2" };
            var speciesData = new Thermo().BuildSpeciesData(speciesList);

            double Tguess = 3800;
            double P = 1000;      
            
            var result = Solve(elements, speciesList, reactants, speciesData, Tguess, P);
            
        }

        public double[][] Solve(string[] elements, string[] species, Dictionary<string, double[]> reactants, Dictionary<string, DataStructure> speciesData, double T, double P)
        {
            double[] Ni = new double[species.Length];
            double[] ni = new double[species.Length];
            double nni = 0;
            double mixtureMassi = 0;

            for (var j = 0; j < species.Length; j++)
            {
                if (reactants.ContainsKey(species[j]))
                {
                    Ni[j] = reactants[species[j]][0];
                }
                else
                {
                    Ni[j] = 0;
                }
                mixtureMassi += Ni[j] * speciesData[species[j]].molwt;
            }

            for (var j = 0; j < species.Length; j++)
            {
                ni[j] = Ni[j] / mixtureMassi;
                nni += ni[j];
            }

            int[][] a = BuildA(elements, species, speciesData);
            double[] bi = BuildB(a, ni, elements, speciesData);

            double[][] propertiesi = new double[species.Length][];
            double hoi = 0;
            double hoiT = 0;
            for (var i = 0; i < species.Length; i++)
            {
                double[] coef;
                if (reactants.ContainsKey(species[i]))
                {
                    if (reactants[species[i]][1] > 1000)
                    {
                        coef = speciesData[species[i]].coefHigh;
                    }
                    else
                    {
                        coef = speciesData[species[i]].coefLow;
                    }
                    propertiesi[i] = Properties(reactants[species[i]][1], P, 1, 1, coef); //would require determining mole fraction of species for multicomponent oxidants/fuels
                    hoi += propertiesi[i][1] * ni[i];
                    hoiT = hoi * reactants[species[i]][1];
                }
            }

            double[] n = new double[species.Length];
            double nn = 0;
            for (var i = 0; i < species.Length; i++)
            {
                
                n[i] = 0.1 / Convert.ToDouble(species.Length);
                nn += n[i];
            }
            //begin iteration calculations
            for (var run = 0; run < 100; run++)
            {
                double[] b = BuildB(a, n, elements, speciesData);

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

                var matrix = BuildMatrix(elements, species, a, b, bi, n, nn, ni, hoiT / T, properties);

                var solution = Gauss(matrix);

                double[] pi = new double[elements.Length];
                Array.Copy(solution, pi, elements.Length);

                double[] delLognj = new double[species.Length];
                double delLogn = solution[elements.Length];
                double delLogT = solution[elements.Length + 1];

                for (var j = 0; j < species.Length; j++)
                {
                    delLognj[j] = 0;
                    for (var i = 0; i < elements.Length; i++)
                    {
                        delLognj[j] += a[i][j] * pi[i];
                    }
                    delLognj[j] += delLogn + properties[j][1] * delLogT - properties[j][3];

                }
                double maxdelLognj = Math.Abs(delLognj[0]);
                for (var i = 1; i < species.Length; i++)
                {
                    if (Math.Abs(delLognj[i]) > maxdelLognj)
                    {
                        maxdelLognj = Math.Abs(delLognj[i]);
                    }
                }

                double lambda1 = 2 / Math.Max(5 * Math.Abs(delLogT), Math.Max(5 * Math.Abs(delLogn), Math.Abs(maxdelLognj)));

                double lambda = Math.Min(1, lambda1);

                for (var j = 0; j < species.Length; j++)
                {
                    n[j] = Math.Exp(Math.Log(n[j]) + lambda * delLognj[j]);
                }
                nn = Math.Exp(Math.Log(nn) + lambda * delLogn);
                T = Math.Exp(Math.Log(T) + lambda * delLogT);

                if (delLogn + maxdelLognj + delLogT < 1e-5)
                {
                    break;
                }


            }
            double[] N = new double[species.Length];
            for (var i = 0; i < species.Length; i++)
            {
                N[i] = n[i] / nn;
            }

            double[][] result = new double[][]
            {
                N,
                new double[] { T },
            };
            return result;
        }

        public double[][] BuildMatrix(string[] elements, string[] species, int[][] a, double[] b, double[] bi, double[] n, double nn, double[] ni, double hoi, double[][] properties)
        {
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
                        matrix[k][elements.Length + 1] += n[j] * properties[j][1] * properties[j][1] + n[j] * properties[j][0];
                    }
                    
                    double ho = 0;
                    for (var i = 0; i < species.Length; i++)
                    {
                        ho += properties[i][1] * n[i];// mixtureMass;
                    }
                    matrix[k][elements.Length + 2] = hoi - ho;
                    //matrix[k][elements.Length + 2] = 0;
                    for (var j = 0; j < species.Length; j++)
                    {
                        matrix[k][elements.Length + 2] += n[j] * properties[j][1] * properties[j][3];
                    }
                }
            }

            return matrix;
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

        public double[] BuildB(int[][] a, double[] n, string[] elements, Dictionary<string, DataStructure> speciesData)
        {
            double[] b = new double[elements.Length];
            for (var i = 0; i < a.Length; i++)
            {
                b[i] = 0;
                for (var j = 0; j < n.Length; j++)
                {
                    b[i] += a[i][j] * n[j];
                }
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
