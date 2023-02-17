
package dke2023;

import java.util.ArrayList;
import java.util.List;

import com.github.chen0040.data.utils.TupleTwo;
import com.github.chen0040.gp.services.RandEngine;
import com.github.chen0040.gp.services.SimpleRandEngine;



public class CollectionUtils {
   public static <T> void shuffle(List<T> a, RandEngine randEngine) {
      for(int i=0; i < a.size(); ++i) {
         int j = randEngine.nextInt(i+1);
         exchange(a, i, j);
      }
   }

   public static <T> void shuffle(List<T> a) {
      RandEngine randEngine = new SimpleRandEngine();
      shuffle(a, randEngine);
   }

   public static <T> void exchange(List<T> a, int i, int j) {
      T temp = a.get(i);
      a.set(i, a.get(j));
      a.set(j, temp);
   }


   public static <T extends Comparable<T>> boolean isBetterThan(T a, T b) {
      return a.compareTo(b) < 0;
   }


   public static <T> TupleTwo<List<T>,List<T>> split(List<T> data, double p) {
      int split_index = (int)(data.size() * p);
      List<T> data1 = new ArrayList<>();
      List<T> data2 = new ArrayList<>();
      for(int i=0; i < split_index; ++i){
         data1.add(data.get(i));
      }
      for(int i=split_index; i < data.size(); ++i) {
         data2.add(data.get(i));
      }
      return new TupleTwo<>(data1, data2);
   }
}