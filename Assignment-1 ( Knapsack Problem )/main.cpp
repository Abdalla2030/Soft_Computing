// Abdalla Fadl Shehata  - 20190305 - 4CS-S3
// Mostafa Abdel Nasser - 20190537 - 4CS-S2
#include <bits/stdc++.h>
using namespace std;
////////////////////////////////////////////////////
#define chrom vector<bool>
int N; // Number of items
int P = 20; // Number of chromosomes   ( Population )
int W; // Size of the knapsack ( Weight of knapsack )
vector<int> weights; // Weights of items
vector<int> values; // Values of items
vector<chrom> population; // Population of chromosomes
int Pc = 70; // Probability of Crossover   (0.4 ~ 0.7)
int Pm = 5; // Probability of Mutation       (0.001 ~ 0.1 )
////////////////////////////////////////////////////
bool randomBoolean(); // Generate random boolean value ( true or false)
void generateInitialPopulation(); // Generates random population of P chromosomes ( 20 Chromosomes )

int generateRandomNumber(int,int);
int calculateFitness(chrom &); // returns the fitness of a chromosome (total values)
pair<int, int> select(); // selects two chromosomes to generate a new generation

pair<chrom, chrom> crossover(pair<int, int> &); // creates the new generation from parents
void mutation(pair<chrom, chrom> &); // mutates the new generation with probability of Pm

bool validChrom(chrom &); // checks the new chromosome is valid or not (total weights <= knapsack weight)
void replacement(pair<chrom, chrom> &); // replaces the new generation with the least fitness chromosomes

pair<int, int> getBestSolution(); // find the best chromosome ( optimal fitness )
bool stop(); // return true when the same result is repeated a thousand times ( same best result keeps occurring )
////////////////////////////////////////////////////
int main(){
    int T;
    cin >> T;
    int counter = 1;
    while (T--) {
        cin >>  W;  // Size of knapsack ( Weight of knapsack )
        cin >> N;  // Number of items
        weights.resize(N);   // Weights Vector
        values.resize(N);      // Values Vector
        for (int i = 0; i < N; i++) {
            cin >> weights[i];
            cin >> values[i];
        }
         generateInitialPopulation();
         while (!stop()) {   // selection & crossover & mutation --> New Population
            auto parents = select();
            auto newGene = crossover(parents);
            mutation(newGene);
            replacement(newGene);
        }
        auto bestSolution = getBestSolution();
        chrom bestChrom = population[bestSolution.second];
        cout<<"===================="<<endl;
         cout<<"==== "<<"Case: #" << counter++ <<" ===="<<endl;
        cout <<" Optimal Value = " << bestSolution.first << endl;
         int numberOFItems = count(bestChrom.begin(), bestChrom.end(), true);
        cout << numberOFItems << endl;
        for (int i = 0; i < N; i++) {
            if (bestChrom[i]){ // if item is selected
                cout << weights[i] <<" " << values[i] << "\n";
            }
        }
    }
}
////////////////////////////////////////////////////
// generate population
bool randomBoolean() { // return 0 or 1 ( true of false )
    return rand() % 2 == 0;
}
void generateInitialPopulation(){
    population.clear();
    while (population.size() < P){
        chrom v(N); // chrom size = number of items
        int totalWeight = 0;
        for (int i = 0; i < N; i++) {
            if (randomBoolean()) { // add this item in knapsack or not
                if (totalWeight + weights[i] <= W) { // make sure this is a feasible solution
                    totalWeight += weights[i];
                    v[i] = true; //  item added in knapsack
                }
            }
        }
        population.push_back(v);
    }
}
////////////////////////////////////////////////////
// fitness function and selection ( Roulette Wheel )
int generateRandomNumber(int mn, int mx){
  return  rand()%(mx-mn + 1) + mn;
}
int calculateFitness(chrom &v) {
    int totalValue = 0;
    for (int i = 0; i < N; i++) {
        totalValue += values[i] * v[i];
    }
    return totalValue;
}
// select 2 old individuals to generate 2 new individuals
// Roulette Wheel
pair<int, int> select(){
    vector<int> cum(P); // Cumulative Fitness
    cum[0] = calculateFitness(population[0]);
    for (int i = 1; i < P; i++) {
        cum[i] = cum[i - 1] + calculateFitness(population[i]);
    }
    int r1, r2;
    int selected1 = -1, selected2 = -1;
    while (selected1 == -1){
        r1 = generateRandomNumber(0, cum.back() - 1); // Generate random number in range of cumulative fitness vector
        if (r1 < cum[0]){ // choose the first chromosome
            selected1 = 0;
            break;
        }
        for (int i = 1; i < P; i++) {
            if (r1 >= cum[i - 1] && r1 < cum[i]){
                selected1 = i; // index of chromosome
                break;
            }
        }
    }
    while (selected2 == -1){
        r2 = generateRandomNumber(0, cum.back() - 1);
        if (r2 < cum[0] && selected1 != 0) { // the second condition to make sure that we not take the same chromosome
            selected2 = 0;
            break;
        }
        for (int i = 1; i < P; i++) {
            if (r2 >= cum[i - 1] && r2 < cum[i] && selected1 != i) {
                selected2 = i; // index of chromosome
                break;
            }
        }
    }
    return make_pair(selected1, selected2);
}
////////////////////////////////////////////////////
// Crossover
pair<chrom, chrom> crossover(pair<int, int> &selected){
    bool Rc = generateRandomNumber(0, 100); // generate random number between 0 and 100
    if (Rc <= Pc){ // Pc = 70
        int Xc = generateRandomNumber(1, N - 1); // generate random number between 1 and number of items -1 to crossover
        chrom p1 = population[selected.first];
        chrom p2 = population[selected.second];
        chrom new1(N), new2(N);
        for (int i = 0; i < Xc; i++){
            new1[i] = p1[i];
            new2[i] = p2[i];
        }
        for (int i = Xc; i < N; i++) {
            new1[i] = p2[i];
            new2[i] = p1[i];
        }
        return make_pair(new1, new2);
    }
    return {};
}
////////////////////////////////////////////////////
// mutation
// flip some bits ( 0->1 and 1->0 )
void mutation(pair<chrom, chrom> &p) {
    for (int i = 0; i < N; i++) {
        if (generateRandomNumber(0, 100) <= Pm){
            p.first[i] = !p.first[i];
        }
    }
    for (int i = 0; i < N; i++) {
        if (generateRandomNumber(0, 100) <= Pm) {
            p.second[i] = !p.second[i];
        }
    }
}
////////////////////////////////////////////////////
// replacement
// Check the new chrom after mutation  satisfy weight constraint or not
bool validChrom(chrom &c) {
    int totalWeight = 0;
    for (int i = 0; i < N; i++) {
        totalWeight += weights[i] * c[i];
    }
    return totalWeight <= W;
}
/* replace 2 chromosomes that were generated after mutation with
the 2 old chromosomes that have minimum fitness values  */
void replacement(pair<chrom, chrom> &p){
    if (validChrom(p.first)){
        int minFitness = 1e8; // initial value
        int minIndex;
        for (int i = 0; i < P; i++) {
            if (calculateFitness(population[i]) < minFitness){
                minFitness = calculateFitness(population[i]);
                minIndex = i;
            }
        }
        if (calculateFitness(p.first) > minFitness) { // replace the old chrom that have minimum fitness value with the new gene
            population[minIndex] = p.first;
        }
    }
    if (validChrom(p.second)) {
        int minFitness = 1e8; // initial value
        int minIndex;
        for (int i = 0; i < P; i++) {
            if (calculateFitness(population[i]) < minFitness) {
                minFitness = calculateFitness(population[i]);
                minIndex = i;
            }
        }
        if (calculateFitness(p.second) > minFitness) {
            population[minIndex] = p.second;// replace the old chrom that have minimum fitness value with the new gene
        }
    }
}
////////////////////////////////////////////////////
// stopping of program
// return the highest fitness value and index of the best chrom
pair<int, int> getBestSolution(){
    int mxFitness = 0;
    int mxIndex = -1;
    for (int i = 0; i < P; i++) {
        if (calculateFitness((population[i])) > mxFitness) {
            mxFitness = calculateFitness(population[i]);
            mxIndex = i;
        }
    }
    return make_pair(mxFitness,mxIndex);
}
// return true when the same result repeated a thousand times
bool stop(){
    static int mxFitness = 0;
    static int occurrence = 0;
    if (getBestSolution().first == mxFitness) {
           occurrence++;
            if (occurrence == 1000){
                    return true;
          }
    }
      else{
        mxFitness = getBestSolution().first;
        occurrence = 1;
    }
    return false;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
