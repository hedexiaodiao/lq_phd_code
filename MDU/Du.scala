import chisel3._
import chisel3.util._

class DivIO(val len: Int) extends Bundle {
  val in = Flipped(DecoupledIO(Vec(2, Output(UInt(len.W)))))
  val flush  = Input(Bool())
  //val sign = Input(Bool())
  val out = ValidIO(Output(UInt((len).W)))
  val rem = Output(UInt((len).W))
  //val res2 = Output(UInt((len).W))
}

class Radix8Divider(len: Int = 64) extends Module {
  val io = IO(new DivIO(len))


  val s_idle :: s_shift :: s_compute ::  Nil = Enum(3)
  val state = RegInit(s_idle)
  //io.in.ready := (state === s_idle)
  val newReq = (state === s_idle) && io.in.valid
  val (a, b) = (io.in.bits(0), io.in.bits(1))
  val shiftReg = Reg(UInt((len * 2).W))
  val (hi,lo) = (shiftReg(len * 2-1, len),shiftReg(len - 1, 0))

  val aReg = RegEnable(a, newReq)
  val bReg = RegEnable(b, newReq)

  val cnt = RegEnable(len.U(log2Up(len+1).W), newReq)//

  when(newReq){
    val canSkipShift = (len.U | Log2(b)) - Log2(a)
    cnt := Mux(canSkipShift >= len.U, len.U, canSkipShift)
    state := s_shift
  }.elsewhen (state === s_shift) {
    shiftReg := aReg << cnt
    state := s_compute
  }.elsewhen (state === s_compute) {
    when(cnt > len.U){
      state := s_idle
    }.elsewhen(cnt === len.U){
      val sr1 = Wire(UInt(len.W))
      val x1enough = hi.asUInt >= bReg.asUInt
      sr1 := Mux(x1enough, hi - bReg, hi)
      shiftReg := Cat(sr1,lo(len-2,0),x1enough)
      state := s_idle
    }.elsewhen(cnt === (len-1).U){
      val sr2 = Wire(UInt(len.W))
      val sr1 = Wire(UInt(len.W))
      val x2enough = hi.asUInt >= bReg.asUInt
      sr2 := Cat(Mux(x2enough, hi - bReg, hi)(len - 2, 0), lo(len-1))
      val x1enough = sr2.asUInt >= bReg.asUInt
      sr1 := Mux(x1enough, sr2 - bReg, sr2)
      shiftReg := Cat(sr1,lo(len-3,0),x2enough,x1enough)
      state := s_idle
    }.elsewhen(cnt === (len-2).U){
      val sr4 = Wire(UInt(len.W))
      val sr2 = Wire(UInt(len.W))
      val sr1 = Wire(UInt(len.W))
      val x4enough = hi.asUInt >= bReg.asUInt
      sr4 := Cat(Mux(x4enough, hi - bReg, hi)(len - 2, 0), lo(len-1))
      val x2enough = sr4.asUInt >= bReg.asUInt
      sr2 := Cat(Mux(x2enough, sr4 - bReg, sr4)(len - 2, 0), lo(len-2))
      val x1enough = sr2.asUInt >= bReg.asUInt
      sr1 := Mux(x1enough, sr2 - bReg, sr2)
      shiftReg := Cat(sr1, lo(len-4,0), x4enough, x2enough, x1enough)
      state := s_idle
    }.otherwise{
      val sr4 = Wire(UInt(len.W))
      val sr2 = Wire(UInt(len.W))
      val sr1 = Wire(UInt(len.W))
      val x4enough = hi.asUInt >= bReg.asUInt
      sr4 := Cat(Mux(x4enough, hi - bReg, hi)(len - 2, 0), lo(len-1))
      val x2enough = sr4.asUInt >= bReg.asUInt
      sr2 := Cat(Mux(x2enough, sr4 - bReg, sr4)(len - 2, 0), lo(len-2))
      val x1enough = sr2.asUInt >= bReg.asUInt
      sr1 := Cat(Mux(x1enough, sr2 - bReg, sr2)(len - 2, 0), lo(len-3))
      shiftReg := Cat(sr1, lo(len-4,0), x4enough, x2enough, x1enough)
      cnt := cnt + 3.U
    }
  }
  val kill = (state=/=s_idle) && io.flush
  when(kill){
    state := s_idle
  }
  io.in.ready := (state === s_idle)
  io.out.valid := RegNext(state) === s_compute && state === s_idle && RegNext(!io.flush)
  io.out.bits := shiftReg(len-1,0)
  //io.quot := shiftReg(len-1,0)
  io.rem := shiftReg(len*2-1,len)//remainder

  //  printf("R8D in.valid %d state %d, a %d b %x\n",io.in.valid,state,a,b)
  //  printf("R8D out.valid %d, out %d rem %d\n",io.out.valid,io.out.bits,io.rem)
  //  printf("R8D cnt %d, in.ready %d\n",cnt,io.in.ready)
  //  printf("????????????????????????????????????????????????????????????????????\n")
}