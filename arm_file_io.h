/*
 * file_io.h
 *
 *  Created on: Dec 7, 2018
 *      Author: gaspard
 */
#include <stdio.h>
#include "xsdps.h" /* SD device driver */
#include "ff.h" /*remember to enable the "xilffs" option in bsp*/
#ifndef SRC_FILE_IO_H_
#define SRC_FILE_IO_H_

volatile u32 FileSize= 800; //Bytes
static FIL fil; /* File object */
static FATFS fatfs;
const char Path[] = "0:/";
UINT NumBytesWritten;
/* Useful links:
 * https://forums.xilinx.com/t5/Embedded-Processor-System-Design/Writing-to-SD-card-using-ff-h/td-p/701398
 * https://forums.xilinx.com/t5/Embedded-Processor-System-Design/ZYNQ-Standalone-BSP-bare-metal-FAT-FS-xilffs-libxilffs-SD-Card/td-p/695332
 * */

/* Note:
 *  - enable the "xilffs" option in bsp
 * 	- increase stack size in "lscript.ld", e.g 0x1000000 represent 16MByte;
 * 	- Only use 8.3 style filename, that is to say 8 characters for name, 3 for extension.
 *
 *
 *
 *
 * */
int readDataset(const char* filename, int* row_ptr_out_size_ptr, int* col_idx_out_size_ptr, int* row_ptr_out, int* col_idx_out){
	FRESULT Res;
	Res = f_mount(&fatfs, Path, 0);
	if (Res!= FR_OK) {
		printf("error to mount\n\r");
		return XST_FAILURE;
	}
	Res = f_open(&fil, filename, FA_READ);
	if (Res!= FR_OK) {
		printf("\rFile  %s  not found\n\r", filename);
		return XST_FAILURE;
	}

	int offset = 0;
	UINT read_bytes;

	Res = f_read(&fil, (void*)row_ptr_out_size_ptr, (UINT)1* sizeof(int),&read_bytes);
	offset += read_bytes;
	xil_printf("\r row_ptr_out_size: bytes to read: %d, real read %d \n\r", 1, read_bytes);
	xil_printf("\r row_ptr_out_size: %d\n\r", *row_ptr_out_size_ptr);
	Res = f_read(&fil, col_idx_out_size_ptr, 1* sizeof(int), &read_bytes);
	offset += read_bytes;
	xil_printf("\r col_idx_out_size: bytes to read: %d, real read %d \n\r", 1, read_bytes);
	xil_printf("\r col_idx_out_size: %d\n\r", *col_idx_out_size_ptr);

	Res = f_read(&fil, row_ptr_out, *row_ptr_out_size_ptr* sizeof(int), &read_bytes);
	offset += read_bytes;
	xil_printf("\r row_ptr_out: bytes to read: %d, real read %d \n\r", *row_ptr_out_size_ptr, read_bytes);

	Res = f_read(&fil, col_idx_out, *col_idx_out_size_ptr* sizeof(int), &read_bytes);
	offset += read_bytes;
	xil_printf("\r col_idx_out: bytes to read: %d, real read %d \n\r", *col_idx_out_size_ptr, read_bytes);
	Res = f_close(&fil);

	xil_printf("total read %d bytes", offset);
	return XST_SUCCESS;
}



int writeDataset(const char* filename, int* row_ptr_out_size_ptr, int* col_idx_out_size_ptr, int* row_ptr_out, int* col_idx_out){
	FRESULT Res;
	Res = f_mount(&fatfs, Path, 0);
	if (Res!= FR_OK) {
		printf("error to mount\n\r");
		return XST_FAILURE;
	}
	Res = f_open(&fil, filename, FA_CREATE_ALWAYS | FA_WRITE);//FA_CREATE_ALWAYS
	if (Res!= FR_OK) {
		printf("\rFile  %s  not found\n\r", filename);
		return XST_FAILURE;
	}

	int offset = 0;
	UINT read_bytes;

	Res = f_write(&fil, (void*)row_ptr_out_size_ptr, (UINT)1* sizeof(int),&read_bytes);
	offset += read_bytes;
	xil_printf("\r row_ptr_out_size: bytes to read: %d, real read %d \n\r", 1, read_bytes);
	xil_printf("\r row_ptr_out_size: %d\n\r", *row_ptr_out_size_ptr);
	Res = f_write(&fil, col_idx_out_size_ptr, 1* sizeof(int), &read_bytes);
	offset += read_bytes;
	xil_printf("\r col_idx_out_size: bytes to read: %d, real read %d \n\r", 1, read_bytes);
	xil_printf("\r col_idx_out_size: %d\n\r", *col_idx_out_size_ptr);

	Res = f_write(&fil, row_ptr_out, *row_ptr_out_size_ptr* sizeof(int), &read_bytes);
	offset += read_bytes;
	xil_printf("\r row_ptr_out: bytes to read: %d, real read %d \n\r", *row_ptr_out_size_ptr, read_bytes);

	Res = f_write(&fil, col_idx_out, *col_idx_out_size_ptr* sizeof(int), &read_bytes);
	offset += read_bytes;
	xil_printf("\r col_idx_out: bytes to read: %d, real read %d \n\r", *col_idx_out_size_ptr, read_bytes);
	Res = f_close(&fil);

	xil_printf("total read %d bytes", offset);
	return XST_SUCCESS;
}

int writeArray(const char* filename, int* array, unsigned int size){
	FRESULT Res;
	Res = f_mount(&fatfs, Path, 0);
	if (Res!= FR_OK) {
		printf("error to mount\n\r");
		return XST_FAILURE;
	}
	Res = f_open(&fil, filename, FA_CREATE_ALWAYS | FA_WRITE);
	if (Res!= FR_OK) {
		printf("\rFile  %s  not found or error to create\n\r", filename);
		return XST_FAILURE;
	}
	Res = f_lseek(&fil, 0); // Pointer to beginning of file
	UINT num_bytes_written =0;
	Res = f_write(&fil, (const void*)array, size* sizeof(int),&num_bytes_written);
	xil_printf("num_bytes_written=%d\r\n", num_bytes_written);
	Res = f_close(&fil);
	if(num_bytes_written == size* sizeof(int))
		return XST_SUCCESS;
	else
		return XST_FAILURE;
}
int readArray(const char* filename, int* array, unsigned int size){
	FRESULT Res;
	Res = f_mount(&fatfs, Path, 0);
	if (Res!= FR_OK) {
		printf("error to mount\n\r");
		return XST_FAILURE;
	}
	Res = f_open(&fil, filename, FA_READ);
	if (Res!= FR_OK) {
		printf("\rFile  %s  not found\n\r", filename);
		return XST_FAILURE;
	}

	UINT num_bytes_read =0;
	Res = f_read(&fil, array, size * sizeof(int),&num_bytes_read);
	xil_printf("num_bytes_read=%d\r\n", num_bytes_read);
	Res = f_close(&fil);
	if(num_bytes_read == size* sizeof(int))
		return XST_SUCCESS;
	else
		return XST_FAILURE;
}



#endif /* SRC_FILE_IO_H_ */