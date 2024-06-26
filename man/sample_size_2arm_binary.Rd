\name{sample_size_2arm_binary}
\alias{sample_size_2arm_binary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the sample size for each arm in a two-arm trial
}
\description{
Calculate the sample size for each arm in a two-arm trial
}
\usage{
sample_size_2arm_binary(p0, p1, p2, diff = 0, selection.prob = 0.9,
                        alpha = 0.1, beta = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p0}{
the successful probability of historical data
  }
  \item{p1}{
the response rate of arm 1
  }
  \item{p2}{
the response rate of arm 2
  }
  \item{diff}{
the equivalence margin
  }
  \item{selection.prob}{
the probability of selection of a superior arm. The default value is selection.prob = 0.9
  }
  \item{alpha}{
the type I error to be controlled. The default value is alpha = 0.1
  }
  \item{beta}{
the type II error to be controlled. The default value is beta = 0.2
  }
}
\value{
\code{sample_size_2arm_binary()} returns required sample size for each arm
}
\references{
Cai, Z., Pan, H., Wu, J., Hsu, C.W. (2024). Uncontrolled Randomized Screening Selection Design for Pediatric Oncology Trials. Accepted in Book Chapter of "Master Protocol Clinical Trial for Efficient Evidence Generation"

Wu, J., Pan, H., & Hsu, C. W. (2022). Two-stage screened selection designs for randomized phase II trials
with time-to-event endpoints. Biometrical Journal, 64(7), 1207-1218

Yap, C., Pettitt, A. & Billingham, L. Screened selection design for randomised phase II oncology trials: an example in chronic lymphocytic leukaemia. BMC Med Res Methodol 13, 87 (2013)
}
\author{
Chia-Wei Hsu, Zongheng Cai, Haitao Pan
}
\examples{
sample_size_2arm_binary(p0 = 0.2, p1 = 0.415, p2 = 0.615)
}
