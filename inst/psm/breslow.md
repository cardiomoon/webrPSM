breslow.dat {robust}	R Documentation

## Breslow Data

## Description

Patients suffering from simple or complex partial seizures were randomized to receive either the antiepileptic drug progabide or a placebo. At each of four successive postrandomization clinic visits, the number of seizures occuring over the previous two weeks was reported.

## Usage

```
data(breslow.dat)
```

## Format

A data frame with 59 observations on the following 12 variables.

`ID`
an integer value specifying the patient identification number.

`Y1`
an integer value, the number of seizures during the first two week period.

`Y2`
an integer value, the number of seizures during the second two week period.

`Y3`
an integer value, the number of seizures during the third two week period.

`Y4`
an integer value, the number of seizures during the fourth two week period.

`Base`
an integer value giving the eight-week baseline seizure count.

`Age`
an integer value giving the age of the parient in years.

`Trt`
the treatment: a factor with levels placebo and progabide.

`Ysum`
an integer value, the sum of Y1, Y2, Y3 and Y4.

`sumY`
an integer value, the sum of Y1, Y2, Y3 and Y4.

`Age10`
a numeric value, Age divided by 10.

`Base4`
a numeric value, Base divided by 4.

## References

Breslow, N. E., and Clayton, D. G. (1993), "Approximate Inference in Generalized Linear Mixed Models," Journal of the American Statistical Association, Vol. 88, No. 421, pp. 9-25.

Thrall, P. F., and Vail, S. C. (1990), "Some Covariance Models for Longitudinal Count Data With Overdispersion," Biometrics, Vol. 46, pp. 657-671.

## Examples

```
data(breslow.dat)
```