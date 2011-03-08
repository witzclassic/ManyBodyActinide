/*
 *  sparseMatrix.h
 *  SparseMatrix
 *
 *  Created by Brad Marston on 4/1/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#import <istream>
#import <string>
#import <cmath>
#import <fstream>

static const double COMPLEX_FRACTION = 0.9;  // may have to increase upwards to 1.0
static const double COMPLEX_CUTOFF = 1.0e-10;

typedef struct {float r, i;} complexEntry; // use floats to save memory
typedef struct {int rowIndex, columnIndex; float value;} realMatrixElement; // use when possible
typedef struct {int rowIndex, columnIndex; complexEntry value;} complexMatrixElement;

class vector;

class sparseMatrix {
    
  friend std::ostream &operator<<(std::ostream &os, sparseMatrix *sm);

    friend long rows(const sparseMatrix *a);
    friend sparseMatrix *readFromMTXFile(const std::string fn);
    friend sparseMatrix *readFromBinaryFile(const std::string fn);
    friend sparseMatrix *sparseMatrixAdd(const sparseMatrix *a, const sparseMatrix *b);
    
    friend void sparseMultiply(const sparseMatrix *h, const vector *v, vector *Av);
    friend void preconditioner(const sparseMatrix *h, double lambda, const vector *q, vector *v);

public:
    
    sparseMatrix(); // the default constructor
    sparseMatrix(const long rows, const long columns, const long realElements, const long complexElements);
    sparseMatrix(const sparseMatrix &rhs);
    ~sparseMatrix();
    bool writeToBinaryFile(const std::string fn);
    void sparseMatrixMultiply(const double);
    
private:
    
    long rows;
    long columns;
    long realElements;
    long complexElements;
    
    realMatrixElement *realOffdiagonalMatrixElement;
    complexMatrixElement *complexOffdiagonalMatrixElement;
    float *diagonalMatrixElement;
    
};

inline sparseMatrix::sparseMatrix()
{
    rows = 0;
    columns = 0;
    realElements = 0;
    complexElements = 0;
    
    realOffdiagonalMatrixElement = NULL;
    complexOffdiagonalMatrixElement = NULL;
    diagonalMatrixElement = NULL;
}

inline sparseMatrix::sparseMatrix(const long _rows, const long _columns, const long _realElements, const long _complexElements)
{
    rows = _rows;
    columns = _columns;
    realElements = _realElements;
    complexElements = _complexElements;

    if (rows != columns) {
        std::cerr << "Not a square matrix\n";
        exit(EXIT_FAILURE);
    }
    
    realOffdiagonalMatrixElement = new realMatrixElement[realElements];
    complexOffdiagonalMatrixElement = new complexMatrixElement[complexElements];
    diagonalMatrixElement = new float[rows];
    
    memset(realOffdiagonalMatrixElement, 0, realElements*sizeof(realMatrixElement));
    memset(complexOffdiagonalMatrixElement, 0, complexElements*sizeof(complexMatrixElement));
    memset(diagonalMatrixElement, 0, rows*sizeof(float));
}

inline sparseMatrix::sparseMatrix(const sparseMatrix &rhs)
{
    rows = (&rhs)->rows;
    columns = (&rhs)->columns;
    complexElements = (&rhs)->complexElements;
    realElements = (&rhs)->realElements;
    
    realOffdiagonalMatrixElement = new realMatrixElement[realElements];
    complexOffdiagonalMatrixElement = new complexMatrixElement[complexElements];    
    diagonalMatrixElement = new float[rows];
    
    memcpy(realOffdiagonalMatrixElement, (&rhs)->realOffdiagonalMatrixElement, realElements*sizeof(realMatrixElement));
    memcpy(complexOffdiagonalMatrixElement, (&rhs)->complexOffdiagonalMatrixElement, complexElements*sizeof(complexMatrixElement));
    memcpy(diagonalMatrixElement, (&rhs)->diagonalMatrixElement, rows*sizeof(float));
}    

inline sparseMatrix::~sparseMatrix()
{
    if (realOffdiagonalMatrixElement) delete[] realOffdiagonalMatrixElement;
    if (complexOffdiagonalMatrixElement) delete[] complexOffdiagonalMatrixElement;
    if (diagonalMatrixElement) delete[] diagonalMatrixElement;
} 
    
inline bool sparseMatrix::writeToBinaryFile(const std::string fn)
{
    std::string binFile = fn + ".bin";
    std::ofstream file(binFile.c_str());
    if (!file) return false;
    
    file.write((char *) &rows, sizeof(long));
    file.write((char *) &columns, sizeof(long));
    file.write((char *) &realElements, sizeof(long));
    file.write((char *) &complexElements, sizeof(long));
    
    file.write((char *) diagonalMatrixElement, rows*sizeof(float));
    
    long offset = realElements * sizeof(realMatrixElement)/4; // hack to work around bug in binary write: break into 4 pieces
    file.write((char *) realOffdiagonalMatrixElement, offset);
    file.write((char *) realOffdiagonalMatrixElement+offset, offset);
    file.write((char *) realOffdiagonalMatrixElement+2*offset, offset);
    file.write((char *) realOffdiagonalMatrixElement+3*offset, offset);    
    
    offset = complexElements * sizeof(complexMatrixElement)/4; // hack to work around bug in binary write: break into 4 pieces
    file.write((char *) complexOffdiagonalMatrixElement, offset);
    file.write((char *) complexOffdiagonalMatrixElement+offset, offset);
    file.write((char *) complexOffdiagonalMatrixElement+2*offset, offset);
    file.write((char *) complexOffdiagonalMatrixElement+3*offset, offset);

    file.close();
    
    return true;
}

inline void sparseMatrix::sparseMatrixMultiply(const double c)
{
    for (long i = 0; i < realElements; i++) {
        realOffdiagonalMatrixElement[i].value = c * realOffdiagonalMatrixElement[i].value;
    }
    
    for (long i = 0; i < complexElements; i++) {
        complexOffdiagonalMatrixElement[i].value.r = c * complexOffdiagonalMatrixElement[i].value.r;
        complexOffdiagonalMatrixElement[i].value.i = c * complexOffdiagonalMatrixElement[i].value.i;
    }
    
    for (long i = 0; i < rows; i++) diagonalMatrixElement[i] = c * diagonalMatrixElement[i];
}


/* friend functions */

inline std::ostream &operator<<(std::ostream &os, sparseMatrix *sm)
{

  os << sm->rows << " " << sm->columns << " " << sm->realElements << " " << sm->complexElements << "\n";
  os << "<============ Diagonal Entries ============>\n";
  for (long i = 0; i < sm->rows; i++) {
      os << i << " " << i << " " <<  sm->diagonalMatrixElement[i] << "\n";
    }
  os << "<============ Real Entries ============>\n";
  for (long i = 0; i < sm->realElements; i++) {
    os << sm->realOffdiagonalMatrixElement[i].rowIndex << " " << sm->realOffdiagonalMatrixElement[i].columnIndex  << " " <<  sm->realOffdiagonalMatrixElement[i].value << "\n";  
    if (sm->realOffdiagonalMatrixElement[i].value == 0)  {
      os << "zero entry! \n";
      break;
    }
  }

  if (sm->complexElements > 0) {
      os << "<============ Complex Entries ============>\n";
      for (long i = 0; i < sm->complexElements; i++) {
	os << sm->complexOffdiagonalMatrixElement[i].rowIndex << " " << sm->complexOffdiagonalMatrixElement[i].columnIndex  << " " <<  sm->complexOffdiagonalMatrixElement[i].value.r << " " << sm->complexOffdiagonalMatrixElement[i].value.i << "\n";
      }
    }
    
  return os;

}

inline long rows(const sparseMatrix *a)
{
    return a->rows;
}
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
inline sparseMatrix *readFromMTXFile(const std::string fn)
{
  
  unsigned long rows, columns, realElements, complexElements;
  unsigned long rowIndex, columnIndex; 
  float realMatrixElement, imagMatrixElement = 0.0;
  unsigned long ir = 0;
  unsigned long ic = 0;

  std::string mtxFile = fn + ".mtx";
  std::ifstream file(mtxFile.c_str());
  if (!file) return NULL;
    
  while (file.peek() == '%') { // strip off comment header
    file.ignore(1024, '\n');
  }
    
  file >> rows >> columns >> realElements >> complexElements;;

  std::cout << fn << " matrix is " << rows << " by " << columns << " with " << realElements << " real " << complexElements << " complex " << std::endl; 
    
  sparseMatrix *sparse = new sparseMatrix(rows, columns, realElements, complexElements); // crease sparse matrix

  //
  // read in matrix elements, sorting into separate real and complex data structures    
  //

  while(file >> rowIndex >> columnIndex >> realMatrixElement >> imagMatrixElement) {
        
    if ((rowIndex < 1) || (rowIndex > rows)) {
      std::cerr << "Row index out of bounds!";
      exit(EXIT_FAILURE);
    }
        
    if ((columnIndex < 1) || (columnIndex > columns)) {
      std::cerr << "Column index out of bounds!";
      exit(EXIT_FAILURE);
    }
        
    rowIndex = rowIndex - 1; // begin with 0, not 1
    columnIndex = columnIndex - 1;
        
    if (rowIndex == columnIndex) 
      sparse->diagonalMatrixElement[rowIndex] += realMatrixElement;

    else {
      if (rowIndex >columnIndex) { // lower half off-diagonal matrix elements are not stored
	std::cerr << "Wrong half for mtx file";
	exit(EXIT_FAILURE);
      }

      if (imagMatrixElement) { // complex entry

	if (ic == complexElements) {
	  std::cout <<  rowIndex << " " << columnIndex << " " << realMatrixElement << " " << imagMatrixElement << std::endl;
	  std::cerr << "Not enough storage allocated for complex matrix elements!";
	  exit(EXIT_FAILURE);
	}            

	sparse->complexOffdiagonalMatrixElement[ic].rowIndex = rowIndex;
	sparse->complexOffdiagonalMatrixElement[ic].columnIndex = columnIndex;
	sparse->complexOffdiagonalMatrixElement[ic].value.r = realMatrixElement;
	sparse->complexOffdiagonalMatrixElement[ic].value.i = imagMatrixElement;
	ic = ic + 1;
      } // complex

      else { // real

	// std::cout <<  ir << " " << ic << " " << rowIndex << " " << columnIndex << " " << realMatrixElement << " " << imagMatrixElement << std::endl;
	if (ir ==  realElements) { // fill in real entries
	  std::cerr << "Not enough storage allocated for real matrix elements!";
	  exit(EXIT_FAILURE);
	}            
            
	sparse->realOffdiagonalMatrixElement[ir].rowIndex = rowIndex;
	sparse->realOffdiagonalMatrixElement[ir].columnIndex = columnIndex;
	sparse->realOffdiagonalMatrixElement[ir].value = realMatrixElement;
	ir = ir + 1;
      } // real
    } // off diag
  } // while real file 
    
  file.close();
    
  //std::cout << sparse;
  return sparse;
}   

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
inline sparseMatrix *readFromBinaryFile(const std::string fn)
{
    std::string binFile = fn + ".bin";
    std::ifstream file(binFile.c_str());
    if (!file) return NULL;
    
    long rows, columns, realElements, complexElements;
    
    file.read((char *) &rows, sizeof(long));
    file.read((char *) &columns, sizeof(long));
    file.read((char *) &realElements, sizeof(long));
    file.read((char *) &complexElements, sizeof(long));


    sparseMatrix *sparse = new sparseMatrix(rows, columns, realElements, complexElements); // crease sparse matrix
    
    file.read((char *) sparse->diagonalMatrixElement, rows*sizeof(float));
    
    long offset = realElements * sizeof(realMatrixElement)/4; // hack to work around bug in binary write: break into 4 pieces
    file.read((char *) sparse->realOffdiagonalMatrixElement, offset);
    file.read((char *) sparse->realOffdiagonalMatrixElement+offset, offset);
    file.read((char *) sparse->realOffdiagonalMatrixElement+2*offset, offset);
    file.read((char *) sparse->realOffdiagonalMatrixElement+3*offset, offset);    
    
    offset = complexElements * sizeof(complexMatrixElement)/4; // hack to work around bug in binary write: break into 4 pieces
    file.read((char *) sparse->complexOffdiagonalMatrixElement, offset);
    file.read((char *) sparse->complexOffdiagonalMatrixElement+offset, offset);
    file.read((char *) sparse->complexOffdiagonalMatrixElement+2*offset, offset);
    file.read((char *) sparse->complexOffdiagonalMatrixElement+3*offset, offset);
    
    file.close();
    
    return sparse;
}

inline sparseMatrix *sparseMatrixAdd(const sparseMatrix *a, const sparseMatrix *b)
{
    if ((a->rows != b->rows) || (a->columns != b->columns)) {
        std::cerr << "Matrices cannot be added because their dimensions differ\n";
        exit(EXIT_FAILURE);
    }
    
    sparseMatrix *c = new sparseMatrix(a->rows, a->columns, a->realElements+b->realElements, a->complexElements+b->complexElements);
    
    long offset = a->realElements;
    memcpy(c->realOffdiagonalMatrixElement, a->realOffdiagonalMatrixElement, a->realElements*sizeof(realMatrixElement));
    memcpy(c->realOffdiagonalMatrixElement+offset, b->realOffdiagonalMatrixElement, b->realElements*sizeof(realMatrixElement));
    
    offset = a->complexElements;
    memcpy(c->complexOffdiagonalMatrixElement, a->complexOffdiagonalMatrixElement, a->complexElements*sizeof(complexMatrixElement));
    memcpy(c->complexOffdiagonalMatrixElement+offset, b->complexOffdiagonalMatrixElement, b->complexElements*sizeof(complexMatrixElement));
        
    for (long i = 0; i < a->rows; i++) {
        c->diagonalMatrixElement[i] = a->diagonalMatrixElement[i] + b->diagonalMatrixElement[i];
    }
    
    return c;
    
}
