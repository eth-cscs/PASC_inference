#ifndef PASC_COMMON_IO_BMP_H
#define	PASC_COMMON_IO_BMP_H

typedef int LONG;
typedef unsigned short WORD;
typedef unsigned int DWORD;

typedef struct tagBITMAPFILEHEADER {
  WORD  bfType;
  DWORD bfSize;
  WORD  bfReserved1;
  WORD  bfReserved2;
  DWORD bfOffBits;
} BITMAPFILEHEADER, *PBITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  DWORD biSize;
  LONG  biWidth;
  LONG  biHeight;
  WORD  biPlanes;
  WORD  biBitCount;
  DWORD biCompression;
  DWORD biSizeImage;
  LONG  biXPelsPerMeter;
  LONG  biYPelsPerMeter;
  DWORD biClrUsed;
  DWORD biClrImportant;
} BITMAPINFOHEADER, *PBITMAPINFOHEADER;

namespace pascinference {

	void bmp_load(std::string filename) {
		std::vector<char> buffer;
		PBITMAPFILEHEADER file_header;
		PBITMAPINFOHEADER info_header;

		std::ifstream file(filename);

		//TODO: control existence
		file.seekg(0,std::ios::end);
		std::streampos length = file.tellg();
		file.seekg(0,std::ios::beg);

		buffer.resize(length);
		file.read(&buffer[0],length);

		file_header = (PBITMAPFILEHEADER)(&buffer[0]);
		info_header = (PBITMAPINFOHEADER)(&buffer[0] + sizeof(BITMAPFILEHEADER));

		coutMaster << buffer[0] << buffer[1] << std::endl;
		coutMaster << file_header->bfSize << std::endl;
		coutMaster << info_header->biWidth << " " << info_header->biHeight << std::endl;
    
	}

	

}

#endif
