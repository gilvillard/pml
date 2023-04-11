#include <iomanip>
#include <NTL/BasicThreadPool.h>
#include <numeric>
#include <random>
#include <NTL/lzz_pX.h>

#include "util.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_interpolant.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* run one bench for specified rdim,cdim,order                */
/*------------------------------------------------------------*/
void one_bench_pmbasis(long rdim, long cdim, long degree)
{
    VecLong shift(degree+1,degree);

    double t1,t2;

 
    double t_pmbasis2x1=-1.0;
    double t_pmbasisgcd = -1.0;
    double t_pmbasisgcd_2 = -1.0;
    double t_pmbasisgcd_3 = -1.0;
            
    // double t_pmbasis_gen_naive=-1.0;
    // double t_pmbasis_gen_middle=-1.0;

    double t_pmbasis_2x1_2=-1.0;
    //double t_pmbasis_generic_middle=-1.0;

    
    double t_gcd_ntl = -1.0;
    zz_pX f0,f1;
    zz_pX p00,p01,p10,p11;
    long s0,s1;
    Mat<zz_pX> pmat;
    pmat.SetDims(2,1);

    Mat<zz_pX> appbas;
    appbas.SetDims(2,2);

    long order;
    cout.precision(3);
    
    if (rdim==2 && cdim==1)
        {
            t_pmbasis2x1=0.0;
            t_pmbasisgcd=0.0;
            random(f0, degree+1);
            random(f1, degree);
            pmat[0][0] = f0;
            pmat[1][0] = f1;
            order = deg(f0) + deg(f1) + 1;
                    
            s0 = shift[0];
            s1 = shift[1];

            // Vincent's
            for (int j = 0; j < 9; j++)
                cout << "\t2x1" << j;
            cout<<endl;
                    
            for (int j = 0; j < 9; j++)
                {
                    t1 = GetWallTime();
                    pmbasis_2x1(p00,p01,p10,p11,f0,f1,order,s0,s1, j, 100);
                    t2 = GetWallTime();
                    t_pmbasis2x1 = t2-t1;
                    cout << "\t" << t_pmbasis2x1;    
                }

            cout<<endl;
            for (int j = 0; j < 9; j++)
                cout << "\tgcd" << j;
            cout << endl;
            // Kevin's
            for (int j = 0; j < 9; j++)
                {
                    t1 = GetWallTime();
                    pmbasis_gcd(appbas, pmat, order, shift, j);
                    t2 = GetWallTime();
                    t_pmbasisgcd = t2-t1;
                    cout << "\t" << t_pmbasisgcd;
                }
            cout<<endl;
            //t1 = GetWallTime();
            //pmbasis_2x1_2(p00,p01,p10,p11,f0,f1,order,s0,s1, 200);
            //t2 = GetWallTime();
            //t_pmbasis_2x1_2 = t2-t1;

            //t1 = GetWallTime();
            //pmbasis_gcd_2(appbas, pmat, order, shift);
            //t2 = GetWallTime();
            //t_pmbasisgcd_2 = t2-t1;
            //t1 = GetWallTime();
            //pmbasis_gcd_3(appbas, pmat, order, shift);
            //t2 = GetWallTime();
            //t_pmbasisgcd_3 = t2-t1;
            //   t1 = GetWallTime();
            // pmbasis_gcd_general_naive(p00,p01,p10,p11,f0,f1,order,s0,s1);
            /// t2 = GetWallTime();
            // t_pmbasis_gen_naive = t2-t1;

            //t1 = GetWallTime();
            //pmbasis_gcd_general_middleprod(p00,p01,p10,p11,f0,f1,order,s0,s1);
            //t2 = GetWallTime();
            //t_pmbasis_gen_middle = t2-t1;

                    
            //t1 = GetWallTime();
            //pmbasis_gcd_generic_middleprod(p00,p01,p10,p11,f0,f1,order,s0,s1);
            //t2 = GetWallTime();
            //t_pmbasis_generic_middle = t2-t1;


            // NTL's
            t1 = GetWallTime();
            GCD(f0, f0, f1);
            t2 = GetWallTime();
            t_gcd_ntl = t2-t1;


            //cout << "\t"<< t_pmbasis_gen_naive << "\t"<< t_pmbasis_gen_middle <<"\t"<< t_pmbasis_generic_middle;
            cout << "gcd ntl" << endl;
            cout << t_gcd_ntl<<endl;
                
        }   
    cout << endl;
}

void multiplication(){
    cout << "degree\t\tnaive\t\tstrassen\t\twaksman\t\tadaptative\t\twaksman2" << endl;

    zz_pX p00, p01, p10, p11, q00,q01,q10,q11, q1, q2, q3, q4, q5, q6, q7, temp;
    VecLong degree_p = {100, 1000, 5000, 11000, 20000, 100000, 200000, 500000,1000000,2000000 };
    double t1,t2, tnaive, tstrassen, twaksman, tadaptative;
    Mat<zz_pX> P,Q;
    zz_pX buf0, buf1, buf2, buf3, buf4;



    for (size_t i = 0; i < degree_p.size(); i++)
        {
            random(p00, degree_p[i]);
            random(p01, degree_p[i]);
            random(p10, degree_p[i]);
            random(p11, degree_p[i]);

            random(q00, degree_p[i]);
            random(q01, degree_p[i]);
            random(q10, degree_p[i]);
            random(q11, degree_p[i]);

            P.SetDims(2,2);
            P[0][0] = p00;
            P[0][1] = p01;
            P[1][0] = p10;
            P[1][1] = p11;
            Q.SetDims(2,2);
            Q[0][0] = q00;
            Q[0][1] = q01;
            Q[1][0] = q10;
            Q[1][1] = q11;
            
            t1 = GetWallTime();
            
            mul(p00, p00, q00);
            mul(temp, p01, q10);
            add(p00, p00, temp);

            mul(temp, p00, q01);
            mul(p01, p01, q11);
            add(p01, p01, temp);

    
            mul(p10, p10, q10);
            mul(temp, p11, q10);
            add(p10, p10, temp);

            mul(temp, p10, q01);
            mul(p11, p11, q11);
            add(p11, p11, temp);

            t2 = GetWallTime();

            tnaive = t2 - t1;

            t1 = GetWallTime();
            add(q1, q00, q10);
            mul(q1, p00, q1);

            add(q2, q01, q11);
            mul(q2, p11, q2);

            sub(q3, p11, p00);
            sub(temp, q10, q01);
            mul(q3, q3, temp);

            sub(q4, p01, p11);
            add(temp, q10, q11);
            mul(q4, q4, temp);

            sub(q5, p01, p00);
            mul(q5, q5, q10);

            sub(q6, p10, p00);
            add(temp, q00, q01);
            mul(q6, q6, temp);

            sub(q7, p10, p11);
            mul(q7, q7, q01);
            
            add(p00, q1, q5);

            add(p01, q2, q3);
            add(p01 , p01, q4);
            sub(p01, p01, q5);

            add(p10, q1, q3);
            add(p10, p10, q6);
            sub(p10, p10, q7);

            add(p11, q2, q7);     
            t2 = GetWallTime();

            tstrassen = t2-t1;
            
            t1 = GetWallTime();
            multiply_waksman(P, P, Q);
            t2 = GetWallTime();
            twaksman = t2-t1;

            t1 = GetWallTime();
            multiply(P, P, Q);
            t2 = GetWallTime();
            tadaptative = t2-t1;

            cout << degree_p[i] << "\t" << tnaive << "\t" << tstrassen << "\t" << twaksman << "\t" << tadaptative << endl ;
        }
    
}
/*------------------------------------------------------------*/
/* run bench on variety of parameters                         */
/*------------------------------------------------------------*/
void run_bench(long nthreads, long nbits, bool fftprime, long rdim=-1, long cdim=-1, long degree=-1, long order=-1)
{
    SetNumThreads(nthreads);

    if (fftprime)
    {
        cout << "Bench pmbasis, FFT prime p = ";
        if (nbits < 25)
        {
            zz_p::UserFFTInit(786433); // 20 bits
            cout << zz_p::modulus() << ", bit length = " << 20 << endl;
        }
        else if (nbits < 35)
        {
            zz_p::UserFFTInit(2013265921); // 31 bits
            cout << zz_p::modulus() << ", bit length = " << 31 << endl;
        }
        else if (nbits < 45)
        {
            zz_p::UserFFTInit(2748779069441); // 42 bits
            cout << zz_p::modulus() << ", bit length = " << 42 << endl;
        }
        else if (nbits < 61)
        {
            zz_p::FFTInit(0);
            std::cout << zz_p::modulus() << ", bit length = " << NumBits(zz_p::modulus()) << std::endl;
        }
        else
        {
            std::cout << "Asking for FFT prime with too large bitsize (> 60). Exiting." << std::endl;
            return;
        }
    }
    else
    {
        cout << "Bench pmbasis, random prime p = ";
        zz_p::init(NTL::GenPrime_long(nbits));
        cout << zz_p::modulus() << ", bit length = " << nbits << endl;
    }

    std::cout << "Note: negative timings for interpolant variants indicate that not enough interpolation points could be found in the base field." << std::endl;
    cout << "rdim\tcdim\tdeg\torder" << endl;
    
    if (rdim==-1) // then cdim==-1 && order==-1, default case
    {
        VecLong szs = {2, 4, 8, 16, 32, 64, 128, 256, 512};

        for (size_t i=0; i<szs.size(); ++i)
        {
            VecLong cdims = {szs[i]/4, szs[i]/2, 3*szs[i]/4};
            for (long j : cdims)
                if (j > 0)
                {
                    long max_order=65536;
                    if (szs[i]==16)
                        max_order=32768;
                    if (szs[i]==32)
                        max_order=16384;
                    if (szs[i]==64)
                        max_order=4096;
                    if (szs[i]==128)
                        max_order=1024;
                    if (szs[i]==256)
                        max_order=256;
                    if (szs[i]==512)
                        max_order=64;
                    for (long k=2; k<=max_order; k=2*k)
                    {
                        one_bench_pmbasis(szs[i],j,k-1); // degree ~ order
                        //one_bench_pmbasis(szs[i],j,k-1,2*k); // degree ~ order/2
                    }
                }
        }
        cout << endl;
    }
    else
        one_bench_pmbasis(rdim,cdim,degree);

}


/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    if (argc!=3 and argc!=7)
        throw std::invalid_argument("Usage: ./time_pmbasis nbits fftprime (rdim cdim degree order)");
    // assume rdim>0 , cdim>0, order>0

    const long nbits = atoi(argv[1]);
    const bool fftprime = (atoi(argv[2])==1);

    if (argc==7)
    {
        const long rdim = atoi(argv[3]);
        const long cdim = atoi(argv[4]);
        const long degree = atoi(argv[5]);
        const long order = atoi(argv[6]);
        run_bench(1,nbits,fftprime,rdim,cdim,degree,order);
    }
    else
        run_bench(1,nbits,fftprime);

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
