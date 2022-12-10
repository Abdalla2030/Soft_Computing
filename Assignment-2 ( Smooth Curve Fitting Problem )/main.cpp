// Abdalla Fadl Shehata  - 20190305 - 4CS-S3
// Mostafa Abdel Nasser - 20190537 - 4CS-S2
#include <bits/stdc++.h>
using namespace std;

#define point pair<double, double>   // point (x,y)
#define chrom vector<double>
#define inverse(x) 1.0 / x

int N; // Number of data points
int D; // Polynomial degree
int P; // Number of chromosomes
int Pc = 70; // Probability of Crossover   (0.4 ~ 0.7)
int Pm = 1; // Probability of Mutation       (0.001 ~ 0.1 )
int t = 0; // current generation
int T = 1000; // Maximum number of generations
int b = 1; // Dependency Factor  ( 1 ~ 5 )

vector<point > points; // the points of the curve
vector<chrom > population; // Population of chromosomes

double randomDouble(double , double );
void generateInitialPopulation(); // Generates random population of P chromosomes

double calculateError(chrom &); // Calculate the  the mean square error
double calculateFitness(chrom &); // returns the fitness of a chromosome (total values)

int generateRandomNumber(int ,int );
pair<chrom, chrom > select(); // selects two chromosomes to generate a new generation

pair<chrom, chrom > crossover(pair<chrom, chrom > &); // creates the new generation from parents

void mutation(pair<chrom, chrom > &); // mutates the new generation with probability of Pm
void mutate(chrom & );

void replacement(pair<chrom, chrom > &); // replaces the new generation with the least fitness chromosomes
pair<double, int> getBestSolution(); // find the best chromosome ( optimal fitness--> minimum error )
////////////////////////////////////////////////////
int main()
{
   ifstream input; // input and output file
    input.open("input.txt");
    ofstream output;
    output.open("output.txt");
    int testCases; // Number of datasets
    input >> testCases;
    int counter = 1;
    while (testCases--){
        input >> N; // Number of data points
        input >> D; // Polynomial Degree
        D++;
        P = max(5,N/2);
        points.clear();
        points.resize(N);
        for (int i = 0; i < N; i++) {
            input >> points[i].first >> points[i].second;
        }
        generateInitialPopulation();
        auto bestSolution  = getBestSolution();
        while (t <= T) {
            auto parents = select();
            auto newG = crossover(parents);
            mutation(newG);
            replacement(newG);
            auto curbestSolution = getBestSolution();
            if (curbestSolution.first > bestSolution.first){
                bestSolution  = curbestSolution;
            }
            t++;
        }
        chrom bestChrom = population[bestSolution.second];
        output<<"===== "<<"Case: #" << counter++ <<" ===== \n";
    output<<" Coefficients :  \n";
        for (double d: bestChrom){
            output << fixed << setprecision(3) << d << " ";
        }
        output << "\nError = " << fixed << setprecision(3) << inverse(bestSolution.first) << "\n";
        output <<"==================== \n";
    }

    return 0;
}
////////////////////////////////////////////////////
// generate population
double randomDouble(double min, double max){ // generate random double number between [-10 : 10 ]
    return min + (double) (rand()) / ((double) (RAND_MAX / (max - min)));
}
void generateInitialPopulation(){
    population.clear();
    while (population.size() < P) {
        chrom v(D); // chrom size =  Polynomial Degree + 1 ( D increased in main function )
        for (int i = 0; i < D; i++) {
            v[i] = randomDouble(-10, 10);
        }
        population.push_back(v);
    }
}
////////////////////////////////////////////////////
// Calculate mean square error and fitness value
double calculateError(chrom &c){   //  calculate mean square error
    double totalSum = 0;
    for (point p: points){
        double x = p.first;
        double y = p.second;
        double sum = 0;
        for (int i = 0; i < D; i++){   // if D = 3 => chromosone size  = 4
            sum += c[i] * pow(x, i); // sum = a0+ a1*x + a2x^2 + a3x^3
        }
        sum -= y; //  - yactual
        sum *= sum; // square
        totalSum += sum;
    }
    totalSum /= N; // 1 / N
    return totalSum;
}
double calculateFitness(chrom &v) {
    return inverse(calculateError(v));
}
////////////////////////////////////////////////////
// selection ( Tournament selection )
int generateRandomNumber(int mn, int mx){
    return mn + (rand() % (mx - mn + 1));
}
pair<chrom, chrom > select() {
    vector<chrom > bestChroms(P);
    for (int i = 0; i < P; i++){
        chrom c1 = population[generateRandomNumber(0, P - 1)];
        chrom c2 = population[generateRandomNumber(0, P - 1)];
        if (calculateFitness(c1) > calculateFitness(c2)) {
            bestChroms[i] = c1;
        } else {
            bestChroms[i] = c2;
        }
    }
    chrom c1 = bestChroms[generateRandomNumber(0, P - 1)];
    chrom c2 = bestChroms[generateRandomNumber(0, P - 1)];
    return make_pair(c1, c2);
}
////////////////////////////////////////////////////
// Crossover ( 2-Point Crossover )
pair<chrom, chrom > crossover(pair<chrom,chrom> &selected){
    int Rc = generateRandomNumber(0, 100); // generate random number between 0 and 100
    if (Rc <= Pc){
        int Xc1 = generateRandomNumber(1,D-1);
        int Xc2 = generateRandomNumber(1,D- 1);
        if (Xc1 > Xc2){ // Xc1 should not equal Xc2
            swap(Xc1, Xc2);
        }
        chrom p1 = selected.first;
        chrom p2 = selected.second;
        chrom new1(D), new2(D);
        for (int i = 0; i < Xc1; i++) {
            new1[i] = p1[i];
            new2[i] = p2[i];
        }
        for (int i = Xc1; i < Xc2; i++){
            new1[i] = p2[i];
            new2[i] = p1[i];
        }
        for (int i = Xc2; i <D; i++) {
            new1[i] = p1[i];
            new2[i] = p2[i];
        }
        return make_pair(new1, new2);
    } else {
        return selected;
    }
}
////////////////////////////////////////////////////
// Mutation ( Non-uniform Mutation )
void mutation(pair<chrom, chrom > &p) {
    mutate(p.first);
    mutate(p.second);
}
void mutate(chrom &c) {
    for (int i = 0; i < D; i++) {
        int Rm = generateRandomNumber(0, 100);
        if (Rm <= Pm){
            double Lx = c[i] - (-10); // delta lower = -10
            double Ux = 10 - c[i];
            double r1 = randomDouble(0, 1);
            double y;
            if (r1 <= 0.5) {
                y = Lx;
            } else {
                y = Ux;
            }
            double r = randomDouble(0, 1);
            double delta = y * (1 - pow(r, (pow(1 - (double) t / T, (double) b))));
            if (y == Lx){
                c[i] -= delta;
            } else {
                c[i] += delta;
            }
        }
    }
}
////////////////////////////////////////////////////
// Replacement ( Elitist Replacement )
void replacement(pair<chrom, chrom > &p){
    int c1 = generateRandomNumber(0, P - 1);
    int c2 = generateRandomNumber(0, P - 1);
    population[c1] = p.first;
    population[c2] = p.second;
}
// return the highest fitness value and index of the best chrom
pair<double, int> getBestSolution(){
    double mx = calculateFitness(population[0]);
    int mxIndex  = 0;
    for (int i = 1; i < P; i++) {
        if (calculateFitness(population[i]) > mx){
            mx = calculateFitness(population[i]);
            mxIndex  = i;
        }
    }
    return make_pair(mx,mxIndex);
}
