\name{optimal_3arm_binary}
\alias{optimal_3arm_binary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find optimal design parameters
}
\description{
Find the optimal parameters used in the \code{get_oc_3arm_binary()} function
}
\usage{
optimal_3arm_binary(p0, p1, p2, p3, alpha = 0.1, beta = 0.2, tot_sample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p0}{
the response rate of historical data
  }
  \item{p1}{
the response rate of arm 1
  }
  \item{p2}{
the response rate of arm 2
  }
  \item{p3}{
the response rate of arm 3
  }
  \item{alpha}{
the type I error to be controlled. The default value is alpha = 0.1
  }
  \item{beta}{
the type II error to be controlled. The default value is beta = 0.2
  }
  \item{tot_sample}{
the required sample size for each arm from function \code{sample_size_3arm_binary()}
  }
}
\value{
\code{optimal_3arm_binary()} returns: (1) alpha: type I error (2) beta: typeII error (3) r1: the maximum number of successes in stage 1 which will terminate trial (4) n1: the number of subjects in stage 1 (5) r2: the maximum number of successes in stage 2 not to warrant further investigation (6) n: the total number of subjects (stage 1 + stage 2) (7) ESS: the expected sample size for each arm (8) PS:the probability of early stopping
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
optimal_3arm_binary(p0 = 0.2, p1 = 0.415, p2 = 0.515, p3 = 0.615,
                    alpha = 0.1, beta = 0.2, tot_sample = 82)
}
