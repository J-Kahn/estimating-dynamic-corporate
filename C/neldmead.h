#include <ctime>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
using namespace std;

double deviation(matrix a, matrix b){
  double dev = 0;
  for(int i = 0; i < a.nrows(); i++){
    dev += abs(a(i) - b(i));
  }
  return dev;
}

void sortby(matrix &mat, matrix &vec){
    imat numers(vec.nrows());
    matrix vec_temp = vec;
    for(int i = 0; i < vec.nrows(); i++){
      numers(i) = i;
    }
    ord ranks = quickOrder(vec_temp, numers, 0, vec.nrows()-1);
    imat order = ranks.index;
    vec_temp = ranks.vect;
    
    //quickSort(vec.point(), 0, vec.nrows() - 1);
    matrix temp = mat;
    for(int i=0;i<vec.nrows();i++){
        vec(i) = vec_temp(i);
        mat.row_sub(temp.row(order(i)), i);
    }
}

class NelderMeadOptimizer {
    public:
      
      double alpha, gamma, rho, sigma;
        matrix vectors;
        matrix values;
        matrix cog;
      int feed, npar, dimension;
        
      NelderMeadOptimizer(int dimension, int npar, double termination_distance=0.001) {
          this->dimension = dimension;
          this->npar = npar;
          srand(time(NULL));
          alpha = 1;
          gamma = 2;
          rho = -0.5;
          sigma = 0.5;
          feed = 0;
          vectors = matrix(dimension + 1, npar);
          cog = matrix(npar);
          values = matrix(dimension + 1);
          this->termination_distance = termination_distance;
      }
      
      void sort(){
        sortby(vectors, values);
      }
      
      // termination criteria: each pair of vectors in the simplex has to
      // have a distance of at most `termination_distance`
      bool done() {
          if (vectors.nrows() < dimension) {
              return false;
          }
          double dev = 100000;
          double maxdev = -100000;
          for (int i=0; i<dimension+1; i++) {
              for (int j=0; j<dimension+1; j++) {
                  if (i==j) continue;
                  dev = deviation(vectors.row(i), vectors.row(j));
                  if(dev > maxdev){
                    maxdev = dev;
                  }
                  if (dev > termination_distance) {
                      cout << maxdev << endl;
                      return false;
                  }
              }
          }
          
          return true;
      }
      

    void get_cog(){
      for(int i = 0; i< npar; i++){
        cog(i) = 0;
      }
      for (int i = 1; i<=dimension; i++) {
        for(int j = 0; j< npar; j++){
          cog(j) = cog(j) + vectors(i, j)  / (double) dimension;
        }
      }
    }
    
    void init(matrix par, double diff){
      vectors.row_sub(par, 0);
      matrix temp = par;
      for (int i = 0; i < dimension; i++){
        temp = par;
        temp(i) = temp(i) + pow(-1.0, i) * diff;
        vectors.row_sub(temp, i+1);
      }
    }
      
      
    matrix reflect(){
      cout << "Relect" << endl;
      //sort(vectors, values);

      
      get_cog();
      

      return cog + (cog - vectors.row(0)) * alpha;
      
    }
        
    matrix step(matrix vec, double score, int &newmax) {
            
      matrix returned;
      
      if(feed < dimension + 1){
        vectors.row_sub(vec, feed);
        values(feed) = score;

        if(feed < dimension){
          feed++;
          returned = vectors.row(feed);
        }
        else{
          feed ++;
          sort();
          returned = reflect();
        }

      }
      else{
        matrix reflected = reflect();

        if(deviation(reflected, vec) == 0.0){
          
          if (score > values(1) && score < values(dimension)) {
              cout << "Intermediate" << endl;
              vectors.row_sub(reflected, 0);
              values(0) = score;
              sort();
              returned = reflect();
              newmax = 0;
          } else if (score > values(dimension)) {
              cout << "Best" << endl;
              vectors.row_sub(reflected, 0);
              values(0) = score;
              returned = cog + (cog - vectors.row(0))*gamma;
              vectors.row_sub(reflected, 0);
              newmax = 1;

            } else {
              cout << "Worst" << endl;
              returned = cog + (cog - vectors.row(0))*rho;
              newmax = 0;
            }
        }
        else{
          if(score > values(0) ){
            cout << "Double step improved" <<endl;
            values(0) = score;
            vectors.row_sub(vec, 0);
            sort();
            reflected = reflect();
            returned = reflected;
          } else{
            if(values(0) > values(1)){
              cout << "Too far" << endl;
                sort();
                reflected = reflect();
                returned = reflected;
                newmax = 0;
            } else {
              cout << "Shrink!!!" << endl;
                feed = 1;
                
                for (int i=0; i<dimension; i++) {
                    vectors.row_sub(vectors.row(dimension) + (vectors.row(i) - vectors.row(dimension))*sigma, i);
                }
                
                matrix temp = vectors.row(0);
                vectors.row_sub(vectors.row(dimension), 0);
                vectors.row_sub(temp, dimension);
                values(0) = values(dimension);
                returned = vectors.row(1);
                
            }
            
          }
        }
      }
      
      return returned;
    }

 
    private:
        double termination_distance;

};
