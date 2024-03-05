# Dani van Enk, 11823526
# primes.py calculates the 1000th prime number
#   it also finds the longest non-prime sequence between primes

# list and variables
primes = [2]
difference = []
non_primes = []
largest_diff = 0
no_of_primes = 1000
nth_prime = 1
nth_number = 3
begin_non_prime = 0
end_non_prime = 0

# Is_prime() returns a prime if it number is a prime
def Is_prime(number):
    # variables
    deviders = 0

    # Check for deviders
    for integer in range(2,number):
        if number%integer == 0:
            deviders += 1
        else:
            pass
    
    # Prime if no deviders
    if deviders == 0:
        return number

# counts all calculated primes and appends them to the primes list
while nth_prime < no_of_primes:
    if Is_prime(nth_number) != None:
        primes.append(Is_prime(nth_number))
        nth_prime += 1

    nth_number += 1

# calculating differences between neighbouring primes
for integer in range(0, len(primes)-2):
    difference.append(primes[integer+1]-primes[integer])

# finding the largest difference and begin and end of non-prime sequence
for integer in range(0, len(difference)-1):
    if largest_diff < difference[integer]:
        largest_diff = difference[integer]
        begin_non_prime = primes[integer] + 1
        end_non_prime = primes[integer+1]

# listing largest non-prime sequence
for integer in range(begin_non_prime, end_non_prime):
    non_primes.append(integer)

# print the result
print("The thousandth prime number is " + str(primes[len(primes)-1]))
print("The next list are the first 1000 primes:")
print(primes)
print("The next list is the longest non-prime sequence under the "
    + str(no_of_primes) + "th prime number, is has a length of "
    + str(len(non_primes)) + ":")
print(non_primes)