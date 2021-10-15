#!/usr/bin/env python3

from struct import unpack
from numpy import (
    array,
    int8,
    uint64,
    uint32,
    uint8,
    zeros,
    )

class FDpacket:
    """
    Encapsulate packet of frequency domain data.
    Header length and structure are based on R2DAQ (aka ArtooDaq)
    """

    BYTES_IN_PAYLOAD = 8192
    BYTES_IN_HEADER = 32
    BYTES_IN_PACKET = BYTES_IN_PAYLOAD + BYTES_IN_HEADER

#    @property
#    def bytes_in_payload(self):
#        return self.bytes_in_payload
#
#    @property
#    def bytes_in_header(self):
#        return self.bytes_in_header

    @property
    def unix_time(self):
        return self._unix_time

    @property
    def pkt_in_batch(self):
        return self._pkt_in_batch

    @property
    def digital_id(self):
        return self._digital_id

    @property
    def if_id(self):
        return self._if_id

    @property
    def user_data_1(self):
        return self._user_data_1

    @property
    def user_data_0(self):
        return self._user_data_0

    @property
    def reserved_0(self):
        return self._reserved_0

    @property
    def reserved_1(self):
        return self._reserved_1

    @property
    def freq_not_time(self):
        return self._freq_not_time

    @property
    def data(self):
        return self._data

    def __init__(self,ut=0,pktnum=0,did=0,ifid=0,ud0=0,ud1=0,res0=0,
    res1=0,fnt=False,data=None):
        """
        Initialize Packet with the given attributes
        """
        # assign attributes
        self._unix_time = ut
        self._pkt_in_batch = pktnum
        self._digital_id = did
        self._if_id = ifid
        self._user_data_0 = ud0
        self._user_data_1 = ud1
        self._reserved_0 = res0
        self._reserved_1 = res1
        self._freq_not_time = fnt
        self._data = data

    def interpret_data(self):
        """
        Returns
        -------
        x : ndarray
            real-valued array represented by the data.
        """
        x = array(self.data, dtype = uint8)
        return x

    @classmethod
    def from_byte_string(cls,bytestr):
        """
        Parse header and data from the given byte string,
        return an object of type Packet
        """
        # check correct size packet
        len_bytes = len(bytestr)
        if not len_bytes == cls.BYTES_IN_PACKET:
            raise ValueError(
            "Packet length is {0} bytes, should have {1} bytes".format(
            len_bytes,cls.BYTES_IN_PACKET))

        # unpack header
        hdr = unpack(">{0}Q".format(cls.BYTES_IN_HEADER//8),
        bytestr[:cls.BYTES_IN_HEADER])

        ut = uint32(hdr[0] & 0xFFFFFFFF)
        pktnum = uint32((hdr[0]>>uint32(32)) & 0xFFFFF)
        #print "Packet number = {0}".format(pktnum)
        did = uint8(hdr[0]>>uint32(52) & 0x3F)
        ifid = uint8(hdr[0]>>uint32(58) & 0x3F)
        ud1 = uint32(hdr[1] & 0xFFFFFFFF)
        ud0 = uint32((hdr[1]>>uint32(32)) & 0xFFFFFFFF)
        res0 = uint64(hdr[2])
        res1 = uint64(hdr[3]&0x7FFFFFFFFFFFFFFF)
        fnt = not (hdr[3]&0x8000000000000000 == 0)
        # print("printing header word 4: hdr[3]", hex(uint64(hdr[3])))
        # pre-allocate data array long enough to hold 8192-byte payload as
        #8-bit integer numbers
        data = zeros(cls.BYTES_IN_PAYLOAD,dtype=uint8)

        # use a holder array of 64-bit data because struct.unpack requires
        #the payload to be unpacked in 64-bit words
        data_64bit = array(unpack(">{0}Q".format(cls.BYTES_IN_PAYLOAD//8),
        bytestr[cls.BYTES_IN_HEADER:]),dtype=uint64)

        #fill data array word by word from data_64bit
        for i in range(len(data_64bit)):
            for j in range(8):
                data[i*8+j] = uint8((data_64bit[i]>>uint64(8*j))&uint64(0xFF))

        return FDpacket(ut,pktnum,did,ifid,ud0,ud1,res0,res1,fnt,data)
