# Statistics III   Midterm Exam # 1    May 19th, 2021

# Use the this notebook as a template to generate a random sample, 
# make the required calculations using the appropriate `Julia` code, 
# and generate a text file with results to be sumitted.

# Type your 9-digit UNAM account number, name, and birth date, as in the following example:
unamID = 123456789
name = "Ana Lucía Ramírez Pedroza"
birth = "26032001"; # ddmmyyyy

# Generate a random sample with the following code:
using Random
n = sum(digits(unamID)[1:3])
Random.seed!(unamID)
sample = randn(n);

include("Stats3.jl")

# Write and/or reuse `Julia` code to obtain an unbiased and consistent 
# point estimation of the probability that a future observation will
# belong to the interval ]0, 1] and save it with the name: point_estimate
point_estimate = Tn("]0,1]", sample)

# Write and/or reuse `Julia` code to obtain a 99% bayesian interval estimation
# of the probability that a future observation will belong to the interval ]0, 1]
# and save it with the name: interval_estimate
r = Bn("]0,1]", sample, 0.99)
interval_estimate = r.intervalo

# Write and/or reuse `Julia` code to obtain a point estimation
# of the 0.6 quantile and save it with the name: xi
xi = cuantil_puntual(0.6, sample)

# Generate a text file using the following code:
filename = string(unamID) * ".txt"
f = open(filename, "w")
info = string.([name, birth, unamID, point_estimate,
       interval_estimate[1], interval_estimate[2], xi, time()])
for i ∈ 1:length(info)
    write(f, info[i])
    write(f, "\n")
end
close(f)

# Submit the text file in Google Classroom
