\name{get_oc_3arm_binary}
\alias{get_oc_3arm_binary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate operating characteristics for Two-Stage Screened Selection Design for Randomized Phase II Trials with Binary Endpoints for 3 arms
}
\description{
Obtain the operating characteristics of Two-Stage Screened Selection Design for Randomized Phase II Trials with Binary Endpoints for 3 arms. The arguments for this function are from outputs of the functions of \code{sample_size_3arm_binary()} and \code{optimal_3arm_binary()}
}
\usage{
get_oc_3arm_binary(r1, r2, n1, n, p0, p1, p2, p3, diff = 0, nsim, seed = 2483)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r1}{
the maximum number of successes in stage 1 which will terminate trial
  }
  \item{r2}{
the maximum number of successes in stage 2 not to warrant further investigation
  }
  \item{n1}{
the number of subjects in stage 1
  }
  \item{n}{
the total number of subjects (stage 1 + stage 2)
  }
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
  \item{diff}{
the equivalence margin. The default value is diff = 0
  }
  \item{nsim}{
the number of simulated trials
  }
  \item{seed}{
the seed. The default value is seed = 2483
  }
}
\value{
\code{get_oc_3arm_binary()} returns: (1) n: total sample size for each arm (2) SSD.Arm.A: selection probability of Arm A (3) SSD.Arm.B: selection probability of Arm B (4) SSD.Arm.C: selection probability of Arm C (5) SSD.No.Arm: the probability of no arms selected (6) diff: the equivalence margin (7) Mean.N.Arm.A: the average number of patients allocated to Arm A (8) Mean.N.Arm.B: the average number of patients allocated to Arm B (9) Mean.N.Arm.C: the average number of patients allocated to Arm C
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
get_oc_3arm_binary(r1 = 4, r2 = 25, n1 = 15, n = 82,
                   p0 = 0.2, p1 = 0.415, p2 = 0.515,
                   p3 = 0.615, nsim = 100)
}
