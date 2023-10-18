import * as fs from "fs";
import * as path from "path";
import * as zlib from "zlib";
import * as readline from "node:readline";
import { Command } from "commander";

import {
  packDiffExPdu,
  unpackDiffExPdu,
  DiffExMode,
} from "../../client/src/util/diffexpdu";

function arrayEq(a1: Uint32Array, a2: Uint32Array): boolean {
  const nElem = a1.length;
  if (a2.length !== nElem) {
    console.log("arrayEq: different lengths");
    return false;
  }
  if (a1.constructor !== a2.constructor) {
    console.log("arrayEq: different constructors");
    return false;
  }
  for (let i = 0; i < nElem; i += 1) {
    if (a1[i] !== a2[i]) {
      console.log("arrayEq: different at position=", i, a1[i], a2[i]);
      return false;
    }
  }
  return true;
}

function drawWithoutReplacement(
  rangeMax: number,
  dst: Uint32Array,
): Uint32Array {
  const alreadyPicked = new Uint8Array(rangeMax);
  for (let i = 0; i < dst.length; i += 1) {
    let choice;
    do {
      choice = (Math.random() * rangeMax) >> 0;
    } while (alreadyPicked[choice]);
    alreadyPicked[choice] = 1;
    dst[i] = choice;
  }
  return dst;
}

function setup(): [Uint32Array, Uint32Array] {
  const n_obs = 3000000;
  const draw = drawWithoutReplacement(n_obs, new Uint32Array(1500000));
  const set1 = draw.subarray(0, ((2 * draw.length) / 3) >> 0);
  const set2 = draw.subarray(((2 * draw.length) / 3) >> 0);
  set1.sort();
  set2.sort();
  return [set1, set2];
}

function warmup(set1: Uint32Array, set2: Uint32Array): void {
  const testCase = {
    mode: DiffExMode.TopN,
    params: { N: 50 },
    set1,
    set2,
  };
  let buf: Uint8Array;
  for (let i = 0; i < 100; i += 1) {
    buf = packDiffExPdu(testCase);
  }
}

function isDirectorySync(path: string): boolean {
  return fs.statSync(path).isDirectory();
}

async function encode(testSpec: string, options: Record<string, string>) {
  if (options.verbose) {
    console.log(`ENCODE, test spec=${testSpec}, bindir=${options.bindir}`);
  }

  // Make output bindir if not already present
  if (!fs.existsSync(options.bindir)) {
    fs.mkdirSync(options.bindir);
  } else if (!isDirectorySync(options.bindir)) {
    console.error("--bindir already exists and is not a directory");
    return;
  }

  const input = testSpec.endsWith(".gz")
    ? fs.createReadStream(testSpec).pipe(zlib.createGunzip())
    : fs.createReadStream(testSpec);
  const rl = readline.createInterface({ input });
  for await (const l of rl) {
    const testCase = JSON.parse(l);
    const buf = packDiffExPdu({
      mode: DiffExMode.TopN,
      params: { N: 50 },
      set1: new Uint32Array(testCase.set1.postings_list),
      set2: new Uint32Array(testCase.set2.postings_list),
    });
    if (options.verbose) console.log(testCase.test_id);
    fs.writeFileSync(path.join(options.bindir, `${testCase.test_id}.bin`), buf);
  }
  rl.close();

  console.log("All encoded.");
}

async function verify(testSpec: string, options: Record<string, string>) {
  if (options.verbose) {
    console.log(`VERIFY, test spec=${testSpec}, binDir=${options.bindir}`);
  }

  if (!isDirectorySync(options.bindir)) {
    console.error("--bindir already exists and is not a directory");
    return;
  }

  const input = testSpec.endsWith(".gz")
    ? fs.createReadStream(testSpec).pipe(zlib.createGunzip())
    : fs.createReadStream(testSpec);
  const rl = readline.createInterface({ input });
  for await (const l of rl) {
    const testCase = JSON.parse(l);
    if (options.verbose) console.log(testCase.test_id);
    const buf = fs.readFileSync(
      path.join(options.bindir, `${testCase.test_id}.bin`),
    );

    const deArgs = unpackDiffExPdu(buf);
    if (deArgs.mode !== DiffExMode.TopN || deArgs.params.N !== 50) {
      console.error(
        `Mimsmatch on mode or params, testid=${testCase.test_id}`,
        testCase,
        deArgs,
      );
      return;
    }
    if (!arrayEq(new Uint32Array(testCase.set1.postings_list), deArgs.set1)) {
      console.error(
        `Set1 mismatch, testid=${testCase.test_id}`,
        testCase,
        deArgs,
      );
      return;
    }
    if (!arrayEq(new Uint32Array(testCase.set2.postings_list), deArgs.set2)) {
      console.error(
        `Set2 mismatch, testid=${testCase.test_id}`,
        testCase,
        deArgs,
      );
      return;
    }
  }

  console.log("All verified.");
}

const program = new Command(); // CLI arg parsing

program
  .command("encode")
  .description("Create block encoded files from test specification.")
  .argument("<test-spec>", "Test specification file name.")
  .option("-v, --verbose", "Verbose output.")
  .requiredOption("-b, --bindir <bindir>", "Encoded binary (output) directory.")
  .action(encode);

program
  .command("verify")
  .description("Validate block encoded files from test specification.")
  .argument("<test-spec>", "Test specification file name.")
  .requiredOption("-b, --bindir <bindir>", "Encoded binary (output) directory.")
  .option("-v, --verbose", "Verbose output.")
  .action(verify);

program.parse();
