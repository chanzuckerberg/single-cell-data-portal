import fs from "fs";

export default function readJson(path: string): any {
  return JSON.parse(fs.readFileSync(path, "utf8"));
}
