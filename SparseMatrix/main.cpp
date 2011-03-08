
#import <fstream>
#import <iostream>

#import "sparseMatrix.h"
#import "vector.h"
#import "davidson.h"

sparseMatrix *readFromFile(std::string fn)
{
    sparseMatrix *sparse;
    
    sparse = readFromBinaryFile(fn);
    if (sparse) { // prefer to read binary file
        std::cout << "Read binary file " << fn + ".bin" << "\n\n";
    } else { // otherwise read text file and write binary file for use next time
        sparse = readFromMTXFile(fn);
        if (sparse) {
            std::cout << "Read MTX file " << fn + ".mtx" << "\n";
            (*sparse).writeToBinaryFile(fn);
            std::cout << "Wrote out binary file " << fn + ".bin" << "\n\n";
        } else {
            std::cerr << fn << " files aren't found.  Exiting...\n";
            exit(EXIT_FAILURE);
        }
    }
    
    return sparse;
}    

const long MAX_ITERATIONS = 100;

long vector::d; 

int main (int argc, char * const argv[]) {
    
    std::cout.precision(8);
    
    std::cout << "\n\nSparse Matrix Diagonalization of Hubbard Hamiltonian\n\n";
    
    srandom(4762); // seed random number generator

    long N = 1;
    std::cout << "Enter number of low-lying eigenvectors desired:  ";
    std::cin >> N;
    double Hubbard_U = 0.0;
    std::cout << "\n\nEnter Hubbard U (eV):   ";
    std::cin >> Hubbard_U;
    std::cout << "U = " << Hubbard_U << "\n\n";
    
    sparseMatrix *U = readFromFile("Coulomb");
    
    sparseMatrix *t = readFromFile("Hopping");
            
    std::cout << "Creating Hamiltonian\n";
    (*U).sparseMatrixMultiply(Hubbard_U);
    sparseMatrix *H = sparseMatrixAdd(U, t);  // H = t + U
    
    std::cout << "Freeing up memory for t and U\n\n";
    delete t;
    delete U;

    std::cout << "Allocating auxilliary memory\n";
    long d = rows(H);
    
    vector::d = d; // set class variable before allocating vectors
    vector *eigenvector = new vector[N];
    vector *v = new vector[DIM];
    vector *Av = new vector[DIM];
    vector *q = new vector;
    
    std::ifstream startVector("LowLyingStates.mtx");
    if (startVector) {
        std::cout << "Reading text file LowLyingStates.mtx\n";
        for (long i = 0; i < N; i++) {
            startVector >> &eigenvector[i];
        }
    } else {
        std::cout << "Creating random start vectors\n";
        for (long i = 0; i < N; i++) {
            random(&eigenvector[i]);
        }
    }
    
    std::cout << "Beginning Davidson-Liu iterations\n";
    long count;
    for (long k = 0; k < N; k++) {
        
        v[0] = eigenvector[k];
        
        bool finished = false;
        count = 0;
        while (!finished && count < MAX_ITERATIONS) {
            finished = davidson(H, v, Av, q, k, eigenvector);
            count += 1;
        }
        
        eigenvector[k] = v[0];
        
    }
    
    delete H;
    delete[] v;
    delete[] Av;
    delete q;
    
    std::cout << count << " iterations. Eigenvalues are:\n";
    for (long i = 0; i < N; i++) {
        std::cout << "\t" << i+1 << "\t" << eigenvalue(&eigenvector[i]) << "\n";
    }
    
    std::cout << "\n\nWriting out low-lying eigenvectors to file LowLyingStates.mtx\n";
    std::ofstream ground("LowLyingStates.mtx");
    ground.precision(10);
    for (long i = 0; i < N; i++) {
        ground << &eigenvector[i];
    }
    ground.flush();
    ground.close();
    delete[] eigenvector;
          
    std::cout << "\n\nFinished!\n";
    std::cout.flush();
    
    return EXIT_SUCCESS;
}
