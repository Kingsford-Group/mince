#ifndef __DECODER_HPP__
#define __DECODER_HPP__

class MinceOpts;

class Decoder {
public:
	void decode(MinceOpts& minceOpts, std::string& ifname, std::string& ofname);
};

#endif //__DECODER_HPP__
