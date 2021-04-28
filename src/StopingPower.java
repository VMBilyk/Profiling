 class NStoppingPower {
/*
nuclear stopping factor for Zr (by Zigler)
 */
   static double ee (double e){
        double a = (32.53*1.008*91.224*(e/1.008))/(1*40*(1.008+91.224)*Math.sqrt(Math.pow(1, (2.0/3.0))+Math.pow(40, (2.0/3.0))));
        return a;

    }
    static double sNN (double e) {
        double a1 = (Math.log(1+ee(e))/(2*(ee(e)+0.10718*Math.pow(NStoppingPower.ee(e), 0.37544))));
        return a1;
    }
    static double sN (double e) {
       double a2 = ((8.462 * 1 * 40 * 1.008 * NStoppingPower.sNN(e))/((1.008 + 91.224) * Math.sqrt(Math.pow(1.0, (2.0/3.0)) + Math.pow(40, (2.0/3.0)))));
       return a2;
    }

}
class EStoppingPower {
    /*
    electrons stopping factor for Zr (by Zigler)
     */
    //Constants from Zigler
    final static double A11=7.603;
    final static double A51=9120;
    final static double A61=405.2;
    final static double A71=0.005765;
    final static double A62=0.02039;
    final static double A72=2704;
    final static double A82=-14.59;
    final static double A92=5.463;
    final static double A102=-0.7333;
    final static double A112=0.04242;
    final static double A122=-0.0008998;

    static double sl(double e) {
        return (A11*Math.pow(e, 0.45));
    }
    static double sh(double e) {
        return ((A51/e)*Math.log(1+A61/e+A71*e));
    }
    static double se1 (double e) {
      return (EStoppingPower.sl(e)*EStoppingPower.sh(e))/(EStoppingPower.sl(e)+EStoppingPower.sh(e));
  }
   /*
   electron stopping power fo energy more then 1000 keV
    */
    //relativistic factor
  static double b (double e) {
      return Math.sqrt(1-1/(Math.pow(((e*0.10626/100000)+1), 2)));
  }

  static double se2(double e) {
      double relF = EStoppingPower.b(e);
      return (A62/(relF*relF) * (Math.log((A72 * relF*relF)/(1 - relF*relF)) - relF*relF - (A82 + A92 * Math.log(e) + A102 * (Math.log(e)*Math.log(e))  + A112 * (Math.log(e)*Math.log(e)*Math.log(e)) + A122 * (Math.log(e)*Math.log(e)*Math.log(e)*Math.log(e)))));
  }
  static double sTp(double e){
      if (e<1000) {
          return se1(e);
      } else {
          return se2(e);
      }
  }
}
class Stopping {
    static double sf(double e) {
        return ((60.22*6.5107/91.224)*(EStoppingPower.sTp(e)+NStoppingPower.sN(e)));
    }
}
class Integral {
    double f (double x){
        return 1/Stopping.sf(x);
    }
    double intRectangle (double x1, double x2, int n) {
        double s, h;
        s = 0;
        h = (x2-x1)/n;
        for (int i = 1; i < n; i++){
            s += ((x1+i*h)-(x1+(i-1)*h))*f(((x1+i*h)+(x1+(i-1)*h))/2);
        }
        return s;

    }
    double intTrapetion(double x1, double x2, int n){
        double s, h;
        s = 0;
        h = (x2-x1)/n;
        for (int i = 1; i < n; i++){
            s += (f(x1+(i-1)*h)+f(x1+i*h))*h/2;
        }
        return s;
    }
    double intSimpson (double x1, double x2, int n){
        int i, z;
        double h,s;

        n=n+n;
        s = f(x1)*f(x2);
        h = (x2-x1)/n;
        z = 4;
        for(i = 1; i<n; i++){
            s = s + z * f(x1+i*h);
            z = 6 - z;
        }
        return (s * h)/3;
    }
}
class Main {
    public static void main(String[] args){
        Integral intS = new Integral();
        //System.out.println(intS.intRectangle(0.0, 1.0,100));
        //System.out.println(intS.intTrapetion(0.0,1.0,100));
        //System.out.println(intS.intSimpson(0,1,1000));
        for (int i = 990; i <= 1000; i++) {
            double x2 = i;

            System.out.println(intS.intTrapetion(1, x2, 1000));
        }


        System.out.println(Stopping.sf(1000));
        //System.out.println(Stopping.sf(999));
    }
}