#include <stdio.h>
#include "platform.h"
#include "xil_printf.h"
#include "main.h"
#include "xtime_l.h"
#include "xparameters.h"



int main() {
	init_platform();
	xil_printf("start. SIZE=%d, F=%d\r\n", SIZE, XPAR_PS7_CORTEXA9_0_CPU_CLK_FREQ_HZ);
	test_battery_arm();
	xil_printf("end\r\n");
	cleanup_platform();
	return 0;
}