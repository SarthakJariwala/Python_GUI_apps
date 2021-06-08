"""
Majority of this PicoHarp Parser was written by skyjur.
Check out the original code here ---
https://github.com/skyjur/picoharp300-curvefit-ui

Modified by Sarthak Jariwala
Fixed CurveHdr (Sarthak)
"""

import datetime
import numpy as np
from ctypes import (
    Structure,
    addressof,
    c_char,
    c_float,
    c_int,
    c_int32,
    c_int64,
    c_uint,
    c_ulong,
    memmove,
    sizeof
)

DISPCURVES = 8
MAXCURVES = 512
MAXCHANNELS = 65536


class tParamStruct(Structure):
    _pack_ = 4
    _fields_ = [
        ("Start", c_float),
        ("Step", c_float),
        ("End", c_float),
    ]


class tCurveMapping(Structure):
    _pack_ = 4
    _fields_ = [
        ("MapTo", c_int),
        ("Show", c_int),
    ]


class TxtHdr(Structure):
    _pack_ = 4
    _fields_ = [
        ("Ident", c_char * 16),
        ("FormatVersion", c_char * 6),
        ("CreatorName", c_char * 18),
        ("CreatorVersion", c_char * 12),
        ("FileTime", c_char * 18),
        ("CRLF", c_char * 2),
        ("CommentField", c_char * 256),
    ]


class BinHdr(Structure):
    _pack_ = 4
    _fields_ = [
        ("Curves", c_int),
        ("BitsPerHistoBin", c_int),
        ("RoutingChannels", c_int),
        ("NumberOfBoards", c_int),
        ("ActiveCurve", c_int),
        ("MeasMode", c_int),
        ("SubMode", c_int),
        ("RangeNo", c_int),
        ("Offset", c_int),
        ("Tacq", c_int),  # in m
        ("StopAt", c_int),
        ("StopOnOvfl", c_int),
        ("Restart", c_int),
        ("DispLinLog", c_int),
        ("DispTimeFrom", c_int),
        ("DispTimeTo", c_int),
        ("DispCountsFrom", c_int),
        ("DispCountsTo", c_int),
        ("DispCurves", tCurveMapping * DISPCURVES),
        ("Params", tParamStruct * 3),
        ("RepeatMode", c_int),
        ("RepeatsPerCurve", c_int),
        ("RepeatTime", c_int),
        ("RepeatWaitTime", c_int),
        ("ScriptName", c_char * 20),
    ]


class BoardHdr(Structure):
    _pack_ = 4
    _fields_ = [
        ("HardwareIdent", c_char * 16),
        ("HardwareVersion", c_char * 8),
        ("HardwareSerial", c_int),
        ("SyncDivider", c_int),
        ("CFDZeroCross0", c_int),
        ("CFDLevel0", c_int),
        ("CFDZeroCross1", c_int),
        ("CFDLevel1", c_int),
        ("Resolution", c_float),
        ("RouterModelCode", c_int),
        ("RouterEnabled", c_int),
        ("RtChan1_InputType;", c_int),
        ("RtChan1_InputLevel", c_int),
        ("RtChan1_InputEdge", c_int),
        ("RtChan1_CFDPresent", c_int),
        ("RtChan1_CFDLevel", c_int),
        ("RtChan1_CFDZeroCross", c_int),
        ("RtChan2_InputType;", c_int),
        ("RtChan2_InputLevel", c_int),
        ("RtChan2_InputEdge", c_int),
        ("RtChan2_CFDPresent", c_int),
        ("RtChan2_CFDLevel", c_int),
        ("RtChan2_CFDZeroCross", c_int),
        ("RtChan3_InputType;", c_int),
        ("RtChan3_InputLevel", c_int),
        ("RtChan3_InputEdge", c_int),
        ("RtChan3_CFDPresent", c_int),
        ("RtChan3_CFDLevel", c_int),
        ("RtChan3_CFDZeroCross", c_int),
        ("RtChan4_InputType;", c_int),
        ("RtChan4_InputLevel", c_int),
        ("RtChan4_InputEdge", c_int),
        ("RtChan4_CFDPresent", c_int),
        ("RtChan4_CFDLevel", c_int),
        ("RtChan4_CFDZeroCross", c_int),
    ]


class CurveHdr(Structure):
    _pack_ = 4
    _fields_ = [
        ("CurveIndex", c_int32),
        ("TimeOfRecording", c_uint),
        ("HardwareIdent", c_char * 16),
        ("HardwareVersion", c_char * 8),
        ("HardwareSerial", c_int32),
        ("SyncDivider", c_int32),
        ("CFDZeroCross0", c_int32),
        ("CFDLevel0", c_int32),
        ("CFDZeroCross1", c_int32),
        ("CFDLevel1", c_int32),
        ("Offset", c_int32),
        ("RoutingChannel", c_int32),
        ("ExtDevices", c_int32),
        ("MeasMode", c_int32),
        ("SubMode", c_int32),
        ("P1", c_float),
        ("P2", c_float),
        ("P3", c_float),
        ("RangeNo", c_int32),
        ("Resolution", c_float),
        ("Channels", c_int32),
        ("Tacq", c_int32),
        ("StopAfter", c_int32),
        ("StopReason", c_int32),
        ("InpRate0", c_int32),
        ("InpRate1", c_int32),
        ("HistCountRate", c_int32),
        ("IntegralCount", c_int64),
        ("reserved", c_int32),
        ("DataOffset", c_int32),
        ("RouterModelCode", c_int32),
        ("RouterEnabled", c_int32),
        ("RtChan_InputType;", c_int32),
        ("RtChan_InputLevel", c_int32),
        ("RtChan_InputEdge", c_int32),
        ("RtChan_CFDPresent", c_int32),
        ("RtChan_CFDLevel", c_int32),
        ("RtChan_CFDZeroCross", c_int32),
    ]


class ParseError(Exception):
    pass


def _read(f, CType):
    data = f.read(sizeof(CType))
    obj = CType()
    memmove(addressof(obj), data, len(data))
    return obj


def _validate_header(header):
    if not header.Ident == "PicoHarp 300" or not header.FormatVersion == "2.0":
        raise ParseError("Does not look like a PicoHarp 300 file.")


class Curve(object):
    res = None
    data = None

    def __repr__(self):
        return "Curve<resolution: %s, size: %s>" % (self.res, len(self.data))


def timefmt(t):
    d = datetime.datetime.fromtimestamp(t)
    return d.strftime("%a %b %d %H:%M:%S %Y")


class PicoharpParser(object):
    _ready = False

    def __init__(self, filename):
        if isinstance(filename, (str)):  # , unicode)):
            filename = open(filename, mode="rb")
        self.f = filename
        self._prepare()

    def _prepare(self):
        self.f.seek(0)

        header = self._header = _read(self.f, TxtHdr)
        """SJ commented this --- it was giving ParseError"""
        #        _validate_header(header)

        bin_header = self._bin_header = _read(self.f, BinHdr)

        self._boards = []
        for i in range(bin_header.NumberOfBoards):
            self._boards.append(_read(self.f, BoardHdr))

        self._curves = []
        for i in range(bin_header.Curves):
            self._curves.append(_read(self.f, CurveHdr))

    def header(self):
        return [(k, getattr(self._header, k)) for k, t in self._header._fields_]

    def no_of_curves(self):
        return self._bin_header.Curves

    def get_curve(self, n):
        header = self._curves[n]
        res = header.Resolution

        self.f.seek(header.DataOffset)
        array = np.fromfile(self.f, c_uint, header.Channels)

        return res, array

    def get_all_curves(self):
        all_curves = []
        for i in range(len(self._curves)):
            all_curves.append(self.get_curve(i))

        return all_curves

    def get_time_window_in_ns(self, curve_no):
        curve = self._curves[curve_no]
        rep_rate = curve.InpRate0
        res, _ = self.get_curve(curve_no)
        time_window_s = (1 / rep_rate) / res  # in seconds

        return time_window_s * 1e9  # in nannoseconds

    def get_integral_counts(self, curve_no):
        curve = self._curves[curve_no]
        integral_counts = curve.IntegralCount

        return integral_counts

    def info(self):
        txthdr = self._header
        binhdr = self._bin_header
        boards = self._boards
        curves = self._curves
        r = []
        w = r.append
        yesno = lambda x: "true" if x else "false"

        w("Ident            : %s" % txthdr.Ident)
        w("Format Version   : %s" % txthdr.FormatVersion)
        w("Creator Name     : %s" % txthdr.CreatorName)
        w("Creator Version  : %s" % txthdr.CreatorVersion)
        w("Time of Creation : %s" % txthdr.FileTime)
        w("File Comment     : %s" % txthdr.CommentField)

        w("No of Curves     : %s" % binhdr.Curves)
        w("Bits per HistoBin: %s" % binhdr.BitsPerHistoBin)
        w("RoutingChannels  : %s" % binhdr.RoutingChannels)
        w("No of Boards     : %s" % binhdr.NumberOfBoards)
        w("Active Curve     : %s" % binhdr.ActiveCurve)
        w("Measurement Mode : %s" % binhdr.MeasMode)
        w("Sub-Mode         : %s" % binhdr.SubMode)
        w("Range No         : %s" % binhdr.RangeNo)
        w("Offset           : %s" % binhdr.Offset)
        w("AcquisitionTime  : %s" % binhdr.Tacq)
        w("Stop at          : %s" % binhdr.StopAt)
        w("Stop on Ovfl.    : %s" % binhdr.StopOnOvfl)
        w("Restart          : %s" % binhdr.Restart)
        w("DispLinLog       : %s" % binhdr.DispLinLog)
        w("DispTimeAxisFrom : %s" % binhdr.DispTimeFrom)
        w("DispTimeAxisTo   : %s" % binhdr.DispTimeTo)
        w("DispCountAxisFrom: %s" % binhdr.DispCountsFrom)
        w("DispCountAxisTo  : %s" % binhdr.DispCountsTo)

        for i in range(DISPCURVES):
            w("---------------------")
            w("Curve No %s" % i)
            w(" MapTo           : %s" % binhdr.DispCurves[i].MapTo)
            w(" Show            : %s" % yesno(binhdr.DispCurves[i].Show))
            w("---------------------")

        for i in range(3):
            w("---------------------")
            w("Parameter No %s" % i)
            w(" Start           : %f" % binhdr.Params[i].Start)
            w(" Step            : %f" % binhdr.Params[i].Step)
            w(" End             : %f" % binhdr.Params[i].End)
            w("---------------------")

        w("Repeat Mode      : %d" % binhdr.RepeatMode)
        w("Repeats per Curve: %d" % binhdr.RepeatsPerCurve)
        w("Repeat Time      : %d" % binhdr.RepeatTime)
        w("Repeat wait Time : %d" % binhdr.RepeatWaitTime)
        w("Script Name      : %s" % binhdr.ScriptName)

        for i, board in enumerate(boards):
            w("---------------------")
            w("Board No %d" % i)
            w(" HardwareIdent   : %s" % board.HardwareIdent)
            w(" HardwareVersion : %s" % board.HardwareVersion)
            w(" HardwareSerial  : %d" % board.HardwareSerial)
            w(" SyncDivider     : %d" % board.SyncDivider)
            w(" CFDZeroCross0   : %d" % board.CFDZeroCross0)
            w(" CFDLevel0       : %d" % board.CFDLevel0)
            w(" CFDZeroCross1   : %d" % board.CFDZeroCross1)
            w(" CFDLevel1       : %d" % board.CFDLevel1)
            w(" Resolution      : %.6f" % board.Resolution)

            if board.RouterModelCode:
                w(" RouterModelCode       : %d" % board.RouterModelCode)
                w(" RouterEnabled         : %d" % board.RouterEnabled)

                w(" RtChan1_InputType     : %d" % board.RtChan1_InputType)
                w(" RtChan1_InputLevel    : %d" % board.RtChan1_InputLevel)
                w(" RtChan1_InputEdge     : %d" % board.RtChan1_InputEdge)
                w(" RtChan1_CFDPresent    : %d" % board.RtChan1_CFDPresent)
                w(" RtChan1_CFDLevel      : %d" % board.RtChan1_CFDLevel)
                w(" RtChan1_CFDZeroCross  : %d" % board.RtChan1_CFDZeroCross)

                w(" RtChan2_InputType     : %d" % board.RtChan2_InputType)
                w(" RtChan2_InputLevel    : %d" % board.RtChan2_InputLevel)
                w(" RtChan2_InputEdge     : %d" % board.RtChan2_InputEdge)
                w(" RtChan2_CFDPresent    : %d" % board.RtChan2_CFDPresent)
                w(" RtChan2_CFDLevel      : %d" % board.RtChan2_CFDLevel)
                w(" RtChan2_CFDZeroCross  : %d" % board.RtChan2_CFDZeroCross)

                w(" RtChan3_InputType     : %d" % board.RtChan3_InputType)
                w(" RtChan3_InputLevel    : %d" % board.RtChan3_InputLevel)
                w(" RtChan3_InputEdge     : %d" % board.RtChan3_InputEdge)
                w(" RtChan3_CFDPresent    : %d" % board.RtChan3_CFDPresent)
                w(" RtChan3_CFDLevel      : %d" % board.RtChan3_CFDLevel)
                w(" RtChan3_CFDZeroCross  : %d" % board.RtChan3_CFDZeroCross)

                w(" RtChan4_InputType     : %d" % board.RtChan4_InputType)
                w(" RtChan4_InputLevel    : %d" % board.RtChan4_InputLevel)
                w(" RtChan4_InputEdge     : %d" % board.RtChan4_InputEdge)
                w(" RtChan4_CFDPresent    : %d" % board.RtChan4_CFDPresent)
                w(" RtChan4_CFDLevel      : %d" % board.RtChan4_CFDLevel)
                w(" RtChan4_CFDZeroCross  : %d" % board.RtChan4_CFDZeroCross)

            w("---------------------")

            for i, curve in enumerate(curves):
                w("---------------------")
                w("Curve Index       : %d" % curve.CurveIndex)
                w("Time of Recording : %s" % timefmt(curve.TimeOfRecording))
                w("HardwareIdent     : %s" % curve.HardwareIdent)
                w("HardwareVersion   : %s" % curve.HardwareVersion)
                w("HardwareSerial    : %d" % curve.HardwareSerial)
                w("SyncDivider       : %d" % curve.SyncDivider)
                w("CFDZeroCross0     : %d" % curve.CFDZeroCross0)
                w("CFDLevel0         : %d" % curve.CFDLevel0)
                w("CFDZeroCross1     : %d" % curve.CFDZeroCross1)
                w("CFDLevel1         : %d" % curve.CFDLevel1)
                w("Offset            : %d" % curve.Offset)
                w("RoutingChannel    : %d" % curve.RoutingChannel)
                w("ExtDevices        : %d" % curve.ExtDevices)
                w("Meas. Mode        : %d" % curve.MeasMode)
                w("Sub-Mode          : %d" % curve.SubMode)
                w("Par. 1            : %f" % curve.P1)
                w("Par. 2            : %.6f" % curve.P2)
                w("Par. 3            : %.6f" % curve.P3)
                w("Range No          : %d" % curve.RangeNo)
                w("Resolution        : %f" % curve.Resolution)
                w("Channels          : %d" % curve.Channels)
                w("Acq. Time         : %d" % curve.Tacq)
                w("Stop after        : %d" % curve.StopAfter)
                w("Stop Reason       : %d" % curve.StopReason)
                w("InpRate0          : %d" % curve.InpRate0)
                w("InpRate1          : %d" % curve.InpRate1)
                w("HistCountRate     : %d" % curve.HistCountRate)
                w("IntegralCount     : %d" % curve.IntegralCount)
                w("reserved          : %d" % curve.reserved)
                w("dataoffset        : %d" % curve.DataOffset)

                if curve.RouterModelCode:
                    w("RouterModelCode      : %d" % curve.RouterModelCode)
                    w("RouterEnabled        : %d" % curve.RouterEnabled)
                    w("RtChan_InputType     : %d" % curve.RtChan_InputType)
                    w("RtChan_InputLevel    : %d" % curve.RtChan_InputLevel)
                    w("RtChan_InputEdge     : %d" % curve.RtChan_InputEdge)
                    w("RtChan_CFDPresent    : %d" % curve.RtChan_CFDPresent)
                    w("RtChan_CFDLevel      : %d" % curve.RtChan_CFDLevel)
                    w("RtChan_CFDZeroCross  : %d" % curve.RtChan_CFDZeroCross)

        return "\n".join(r)
