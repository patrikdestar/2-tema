#ifndef COMPLEX_H
#define COMPLEX_H

#include <iostream>

/**
 * Compleksiniu skaiciu klase
*/
class Complex {
    private:
        double real, imag;

    public:
        Complex(double r = 0, double i = 0)
        {
            real = r;
            imag = i;
        }
        void print() {
            if (real != 0){
                if (imag > 0){
                    std::cout << real << "+" << imag << "i";
                } else if (imag < 0){
                    std::cout << real << "" << imag << "i";
                } else if (imag == 0){
                    std::cout << real;
                }
            }else {
                if (imag > 0){
                    std::cout << imag << "i";
                } else if (imag < 0){
                    std::cout << imag << "i";
                } else if (imag == 0){
                    std::cout << imag;
                }
            }
            }
        // The global operator function is made friend of this
        // class so that it can access private members
        friend Complex operator+(Complex const& c1, Complex const& c2);
        friend Complex operator-(Complex const& c1, Complex const& c2);
        friend Complex operator*(Complex const& c1, Complex const& c2);
};

Complex operator+(Complex const& c1, Complex const& c2)
{
    return Complex(c1.real + c2.real, c1.imag + c2.imag);
}
Complex operator-(Complex const& c1, Complex const& c2)
{
    return Complex(c1.real - c2.real, c1.imag - c2.imag);
}
Complex operator*(Complex const& c1, Complex const& c2)
{
    // c1.real*c1.real + c1.real*c2.imag + c1.imag*c2.real + c1.imag*c2.imag;
    return Complex(c1.real*c2.real - c1.imag*c2.imag, c1.real*c2.imag + c2.real*c1.imag);
}


#endif