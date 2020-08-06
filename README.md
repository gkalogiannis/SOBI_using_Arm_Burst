# SOBI_using_Arm_Burst
Designing and prototyping a Brain-Computer Interface (BCI) using embedded systems is mainly an
integration process where designers connect a set of custom IP cores together using standard buses in
order to build their custom systems. Such IPs are commonly modelling central functionalities of a BCI
system such as Electroencephalography (EEG) signal capturing, preprocessing, filtering, de-mixing, feature
extraction, classification and mapping.
System-On-Chips (SoC’s), as part of embedded systems, hosts most of the necessary computing
components into a single chip. Writing applications for the latter is focused and dictated by the BCI
application needs. In such cases fast computation must be served in both hardware and software level.
The purpose of this study is to exploit on how to improve the performance of BCI SoC applications in
software level by using Arm Assembly features such the data transfer burst technique. We have chosen
and accelerate the demanding algorithm of Second Order Blind Identification (SOBI) commonly used for
signal separation in BCI systems that are dedicated in recognizing human motor imagery movements. In
our BCI system imagery motor movements are recognized as mu and beta event-related
desynchronization (ERD) and event-related synchronization (ERS) patterns inside EEG recordings. The
algorithm is implemented using C and Arm assembly while measurement of execution time is performed
using Arm Keil MDK, targeting the STM32 Nucleo Cortex-M0 and Cortex-M4 developing boards.
