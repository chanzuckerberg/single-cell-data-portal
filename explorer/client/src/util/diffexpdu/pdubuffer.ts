export class PduBuffer {
  private pos: number; // pointer to unused

  private buf: Uint8Array; // a view into underlying buffer

  private view: DataView; // dataview on the uint8array view

  constructor(capacity = 1024) {
    this.pos = 0;
    this.buf = new Uint8Array(capacity);
    this.view = new DataView(this.buf.buffer);
  }

  ensureFreeCapacity(freeSpaceNeeded: number): void {
    const requiredCapacity = this.pos + freeSpaceNeeded;
    if (requiredCapacity < this.buf.byteLength) return;

    // else, expand
    if (requiredCapacity > 0x80000000) {
      throw new Error("May not grow buffer beyond 2GiB");
    }
    let newCapacity = this.buf.byteLength;
    do {
      newCapacity <<= 1; // eslint-disable-line no-bitwise -- meaningless error!
    } while (requiredCapacity > newCapacity);

    const newBuf = new Uint8Array(newCapacity);
    newBuf.set(this.buf);
    this.buf = newBuf;
    this.view = new DataView(this.buf.buffer);
  }

  advance(offset: number): void {
    this.pos += offset;
  }

  asUint8Array(): Uint8Array {
    return new Uint8Array(this.buf.buffer, this.buf.byteOffset, this.pos);
  }

  asUint32Array(): Uint32Array {
    return new Uint32Array(this.buf.buffer, this.buf.byteOffset, this.pos / 4);
  }

  addUint8(value: number): void {
    this.view.setUint8(this.pos, value);
    this.pos += 1;
  }

  addUint16(value: number): void {
    this.view.setUint16(this.pos, value, true);
    this.pos += 2;
  }

  addUint32(value: number): void {
    this.view.setUint32(this.pos, value, true);
    this.pos += 4;
  }

  add(bytes: Uint8Array): void {
    this.buf.set(bytes, this.pos);
    this.pos += bytes.byteLength;
  }
}
