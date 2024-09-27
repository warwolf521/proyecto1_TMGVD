#include "MurmurHash3.h"
#include <iostream>
#include <string>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <unordered_set>
#include <bitset>
#include <iomanip>

int b = 16;
int m = 1 << b;
double alpha = 0.72129999569229111;

using namespace std;

//funcion para calcular el hash de 64 bits usado en HyperLogLog
uint64_t h(const void* key, int len) {
    uint64_t out[2];
    MurmurHash3_x64_128(key, len, 0, out);
    return out[0];
}

void add(string str, vector<unsigned int>* M){

    uint64_t x = h(str.c_str(), str.length());   //se calcula el hash del k-mer o minimizer

    int index = x >> (sizeof(x) * 8 - b);  //se obtiene los b bits mas significativos para seleccionar la posicion en el sketch

    uint64_t p = x << b;  //se obtienen los 64-b bits para obtener los leading zeros

    (*M)[index] = max( (*M)[index], static_cast<unsigned int>(__builtin_clzll(p) + 1));  //se agrega la nueva observacion al sketch

}

unsigned long estimate(vector<unsigned int>* M){

    double Z = 0;

    for (int i = 0; i < m; i++){    // Calculamos la media armonica del Arreglo M
        Z += 1.0 / (1 << (*M)[i]);   // Z = sumatoria (2^-M[j])
    }

    Z = 1.0 / Z;                    //  Z = Z^-1
    
    unsigned long E = alpha * pow((unsigned long)m, 2)* Z;

    if (E <= (5.0 / 2.0) * m) {
        int V = count((*M).begin(), (*M).end(), 0);  // conteo de registros en cero
        if (V > 0) {
            E = m * log(static_cast<double>(m) / V);  // correccion para rango menor
        }
    }

    return static_cast<unsigned long>(E);
    
}

//funcion que retorna todos los k-mers de una secuencia
vector<string> get_kmers(string sequence, int k){
    vector<string> kmers;
    int size = sequence.size() - k +1;

    for(int i = 0;i<size;i++){
        kmers.push_back(sequence.substr(i,k));
    }

    return kmers;
}

//funcion para obtener el sketch M basado en k-mers de una secuencia dada
vector<unsigned int> k_mers(string input, int k){
    vector<unsigned int> M(1 << b, 0);
    vector<string> kmers = get_kmers(input, k);
    for (auto kmer : kmers) {
        add(kmer, &M);
    }
    
    return M;
}

//Esta funcion para calcular todos los minimizer de sucuencia dado un k y w
vector<string> get_minimizers(string sequence, int w, int k){

    vector<string> min_kmers;
    int windows = sequence.size() - w + 1;
    for(int i = 0;i <windows; i++){
        string minimizer = sequence.substr(i,k);
        for(int j = i+1;j<i+w-k+1;j++){
            
            string s = sequence.substr(j,k);

            if(s.compare(minimizer) < 0){
                minimizer = s;
            }
        }
        min_kmers.push_back(minimizer);
    }
    // iteramos a traves de la palabra, obtieniendo cada k_mers, desde inicio hasta w
    // el limite de i es w-k+1 ya que siempre cortaremos el string desde i hasta i+k
    
    return min_kmers;
}

//funcion para obtener el sketch M basado en minimizers de una secuencia dada
vector<unsigned int> minimizers(string str, int w, int k){
    //definimos la estructura
    vector<unsigned int> M(1 << b, 0);
    vector<string> min = get_minimizers(str, w, k);

    for (auto x : min) {
        add(x, &M);
    }

    return M;
}

//funcion exacta para obtener la cardinalidad de los k-mers distintos de una secuencia
unordered_set<string> exact_kmer(string sequence, int k) {
    unordered_set<string> unique_kmers;
    vector<string> kmers;
    for (size_t i = 0; i <= sequence.length() - k; ++i) {
        kmers.push_back(sequence.substr(i, k));
    }

    for (const auto& kmer : kmers) {
        unique_kmers.insert(kmer);
    }

    return unique_kmers;
}

//funcion exacta para obtener la cardinalidad de los minimizers distintos de una secuencia
unordered_set<string> exact_minimizer(string sequence, int w, int k){
    unordered_set<string> unique_min;
    vector<string> min_kmers = get_minimizers(sequence, w, k);

    for (const auto& kmer : min_kmers) {
        unique_min.insert(kmer);
    }

    return unique_min;
}

//funcion para calcular el ERM y EAM dado dos vectores que contengan los jaccard estimados y jaccard exactos
void calcError(vector<float> J_exact, vector<float> J_hll){

    int n = J_exact.size(), zeros = 0;
    if(J_exact.size() != J_hll.size()){
        cout << "error cardinalidades" << endl;
        return;
    }
    float ERM = 0, EAM = 0;

    for (int i = 0; i < n; i++){
        EAM += abs(J_hll[i] - J_exact[i]);
        if(J_exact[i] == 0.0){
            zeros++;
            continue;
        }
        ERM += (abs(J_hll[i] - J_exact[i]) / J_exact[i]);
    }

    cout << "Error Relativo Medio: " << ERM / (float)(n-zeros) << " Error Absoluto Medio: " << EAM / (float)n << endl;
}

//funcion para calcular jaccard estimado entre 2 secuencias usando k-mers
float jaccard_kmers(string A, string B, int k){
    vector<unsigned int> M_A =  k_mers(A, k);
    vector<unsigned int> M_B =  k_mers(B, k);
    vector<unsigned int> M_AUB((1 << b), 0);

    for (int i = 0; i < (1 << b); i++){
        M_AUB[i] = max(M_A[i], M_B[i]);
    }


    unsigned long est_A = estimate(&M_A);
    unsigned long est_B = estimate(&M_B);
    unsigned long est_AUB = estimate(&M_AUB);
    
    return ((float)est_A + (float)est_B - (float)est_AUB) / (float)est_AUB;

}

//funcion para calcular jaccard estimado entre 2 secuencias usando minimizers
float jaccard_minimizers(string A, string B, int w, int k){

    vector<unsigned int> M_A =  minimizers(A, w, k);
    vector<unsigned int> M_B =  minimizers(B, w, k);
    vector<unsigned int> M_AUB((1 << b), 0);

    for (int i = 0; i < (1 << b); i++){
        M_AUB[i] = max(M_A[i], M_B[i]);
    }
    unsigned long est_A = estimate(&M_A);
    unsigned long est_B = estimate(&M_B);
    unsigned long est_AUB = estimate(&M_AUB);

    return ((float)est_A + (float)est_B - (float)est_AUB) / (float)est_AUB;

}

//funcion para calcular jaccard exacto entre 2 secuencias usando k-mers
float exact_jaccard_kmers(string A, string B, int k){
    unordered_set<string> M_A =  exact_kmer(A, k);
    unordered_set<string> M_B =  exact_kmer(B, k);
    unordered_set<string> M_AUB;

    for(auto x : M_A){
        M_AUB.insert(x);
    }
    for(auto x : M_B){
        M_AUB.insert(x);
    }

    unsigned long est_A = M_A.size();
    unsigned long est_B = M_B.size();
    unsigned long est_AUB = M_AUB.size();
    
    return ((float)est_A + (float)est_B - (float)est_AUB) / (float)est_AUB;

}

//funcion para calcular jaccard exacto entre 2 secuencias usando minimizers
float exact_jaccard_minimizers(string A, string B, int w, int k){

    unordered_set<string> M_A =  exact_minimizer(A, w, k);
    unordered_set<string> M_B =  exact_minimizer(B, w, k);
    unordered_set<string> M_AUB;

    for(auto x : M_A){
        M_AUB.insert(x);
    }
    for(auto x : M_B){
        M_AUB.insert(x);
    }

    unsigned long est_A = M_A.size();
    unsigned long est_B = M_B.size();
    unsigned long est_AUB = M_AUB.size();
    
    return ((float)est_A + (float)est_B - (float)est_AUB) / (float)est_AUB;
}


int main() {
    
    bool calculo_exacto = true; // Cálculos exactos son lentos, esta variable permite ver algunos resultados con un tiempo de ejecución menor
                                 // False -> cálculo de cardinalidades y similitud jaccard hyperloglog.
                                 // True -> Lo anterior + cálculo de cardinalidades exactas, similitud jaccard exacta 
                                 // y cálculo de error entre jaccard hyperloglog y exacto.
                                 
    string str;
    vector<string> dna_sequences(5,"");
    vector<int> Ks = {20, 25};
    vector<int> Ws= {50, 100};
    
    ifstream myfile("GCF_000373685.1_ASM37368v1_genomic.fna");
    getline(myfile, str);
    
    int index = 0;
    while(getline(myfile, str)){
        for(int i=0;i<str.size();i++){
            if(str[i] == 'N')str.erase(i,1);
        }
        if(str[0] == '>'){
            index++;
            if(index == 5)break;
            continue;
        }
        dna_sequences[index].append(str);
    }

    vector<float> vector_jaccard_kmers;
    vector<float> vector_exact_jaccard_kmers;

    vector<float> vector_jaccard_minimizers;
    vector<float> vector_exact_jaccard_minimizers;

    cout << "cardinalidades k-mers" << endl;
    if(calculo_exacto){
        for(auto k : Ks){
            cout << "k: " << k << endl;
            for(int i=0;i<dna_sequences.size();i++){
                vector<unsigned int> estKmers = k_mers(dna_sequences[i], k);
                cout << "sequence: " << i << "     estimada: " << estimate(&estKmers) << " exacta: " << exact_kmer(dna_sequences[i], k).size() << endl;
            }
        }
    }else{
        for(auto k : Ks){
            cout << "k: " << k << endl;
            for(int i=0;i<dna_sequences.size();i++){
                vector<unsigned int> estKmers = k_mers(dna_sequences[i], k);
                cout << "sequence: " << i << " estimada: " << estimate(&estKmers) << endl;
            }
        }
    }
    cout << endl << "jaccard k-mers" << endl;
    if(calculo_exacto){
        for(auto k : Ks){
            cout << "k: " << k << endl;
            for(int i=0;i<dna_sequences.size();i++){
                for(int j=i+1;j<dna_sequences.size();j++){
                    
                    float jaccardKmer = jaccard_kmers(dna_sequences[i], dna_sequences[j], k);
                    float exactJaccardKmer = exact_jaccard_kmers(dna_sequences[i], dna_sequences[j], k);

                    vector_jaccard_kmers.push_back(jaccardKmer);
                    vector_exact_jaccard_kmers.push_back(exactJaccardKmer);
                    
                    cout << setprecision(6) << fixed << "sequence " << i << " and " << j << "     estimada: " << jaccardKmer  << " exacta: " << exactJaccardKmer << endl;
                }
            }
        }
        calcError(vector_exact_jaccard_kmers, vector_jaccard_kmers);
    }else{
        for(auto k : Ks){
            cout << "k: " << k << endl;
            for(int i=0;i<dna_sequences.size();i++){
                for(int j=i+1;j<dna_sequences.size();j++){
                    cout << setprecision(6) << fixed << "sequence " << i << " and " << j << " estimada: " << jaccard_kmers(dna_sequences[i], dna_sequences[j], k) << endl;
                }
            }
        }
    }
    

    cout << endl << "cardinalidades minimizers" << endl;
    if(calculo_exacto){
        for(auto k : Ks){
            for(auto w : Ws){
                cout << "k: " << k << " w: " << w << endl;
                for(int i=0;i<dna_sequences.size();i++){
                    vector<unsigned int> estMinimizer = minimizers(dna_sequences[i], w, k);
                    cout << "sequence: " << i << "     estimada: " << estimate(&estMinimizer) << " exacta: " << exact_minimizer(dna_sequences[i], w, k).size() << endl;
                }
            }
        }
    }else{
        for(auto k : Ks){
            for(auto w : Ws){
                cout << "k: " << k << " w: " << w << endl;
                for(int i=0;i<dna_sequences.size();i++){
                    vector<unsigned int> estMinimizer = minimizers(dna_sequences[i], w, k);
                    cout << "sequence: " << i << " estimada: " << estimate(&estMinimizer) << endl;
                }
            }
        }
    }

    cout << endl << "jaccard minimizers" << endl;
    if(calculo_exacto){
        for(auto k : Ks){
            for(auto w : Ws){
                cout << "k: " << k << " w: " << w << endl;
                for(int i=0;i<dna_sequences.size();i++){
                    for(int j=i+1;j<dna_sequences.size();j++){

                        float jaccardMinimizer = jaccard_minimizers(dna_sequences[i], dna_sequences[j], w, k);
                        float exactJaccardMinimizer = exact_jaccard_minimizers(dna_sequences[i], dna_sequences[j], w, k);

                        vector_jaccard_minimizers.push_back(jaccardMinimizer);
                        vector_exact_jaccard_minimizers.push_back(exactJaccardMinimizer);
                        
                        cout << setprecision(6) << fixed << "sequence " << i << " and " << j << "     estimada: " << jaccardMinimizer  << " exacta: " << exactJaccardMinimizer << endl;
                    }
                }
            }
            
        }
        calcError(vector_exact_jaccard_minimizers, vector_jaccard_minimizers);
    }else{
        for(auto k : Ks){
            for(auto w : Ws){
                cout << "k: " << k << " w: " << w << endl;
                for(int i=0;i<dna_sequences.size();i++){
                    for(int j=i+1;j<dna_sequences.size();j++){
                        cout << setprecision(6) << fixed << "sequence " << i << " and " << j << " estimada: " << jaccard_minimizers(dna_sequences[i], dna_sequences[j], w, k) << endl;
                    }
                }
            }
            
        }
    }
    
    

    return 0;
}
