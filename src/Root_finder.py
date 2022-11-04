#! python3

__author__ = "Simanta Barman"
__email__ = "barma017@umn.edu"

import logging

class Root_finder:
    def __init__(self, function):
        self.function = function

    def bisect(self, a=0, b=1, precision=0.00001, iteration=10):
        """Returns the midpoint of the inerval where the optimal solution lies."""
        logging.info(f'Using Bisection Method:'
                     f"\n{'Iteration':>9}{'a':>16}{'b':>16}{'midpoint':>16}{'dk':>16}")
        # Start of Bisection Algorithm
        # Step 0: Initializing, k, a0, b0
        k, a, b = 1, a, b
        while True:
            # Step 1: Find derivative at midpoint, dk = f'(a+b/2)
            midpoint = (a+b) / 2
            # Inputting differentiated function so, not performing differentiation again here.
            dk = float(self.function(midpoint))
            logging.info(f'{k:>9}{a:>16n}{b:>16n}{midpoint:>16n}{dk:>16n}')
            # Step 2:
            a, b = (a, midpoint) if dk > 0 else (midpoint, b)
            # Step 3: Terminating condition check and terminate or return to 1
            if abs(b - a) < precision or k == iteration:
                logging.info('')
                return (a+b)/2
            k += 1
    
    def newton_raphson(self, x=1, precision=0.00001, iteration=8):
        """Returns the root using Newton-Raphson Method"""
        logging.info(f'Using Newton-Raphson Method:'
                f"\n{'Iteration':>9}{'xᵢ':>32}{'f(xᵢ)':>32}{'f_prime(xᵢ)':>32}{'xᵢ₊₁':>32}")
        
        # Start of Bisection Algorithm
        # Step 0: Initializing k, xᵢ
        k, xᵢ = 1, x

        def ddx(func, x):
            h = 1e-6
            return (func(x+h) - func(x-h)) / (2*h)
        
        while True:
            # Step 1: Find the f(xᵢ)
            f_xᵢ = self.function(xᵢ)
            # Step 2: Find the derivative, f'(xᵢ)
            f_prime_xᵢ = ddx(self.function, xᵢ)
            # Step 3: Find the next, xᵢ₊₁
            xᵢ_new = xᵢ -  f_xᵢ / f_prime_xᵢ
            logging.info(f"{k:>9}{xᵢ:>32}{f_xᵢ:>32}{f_prime_xᵢ:>32}{xᵢ_new:>32}")
            # Step 4: Terminating condition check and terminate or return to 1
            if abs(xᵢ_new - xᵢ) < precision or k == iteration:
                logging.info('')
                return xᵢ_new
            
            k += 1
            xᵢ = xᵢ_new
