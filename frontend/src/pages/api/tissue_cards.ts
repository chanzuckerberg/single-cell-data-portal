import type { NextApiRequest, NextApiResponse } from "next";
import readJson from "src/common/utils/api/readJson";

const data = readJson("src/views/CellCards/common/fixtures/allTissues.json");

async function handler(_: NextApiRequest, res: NextApiResponse) {
  try {
    res.status(200).json(data);
  } catch (error) {
    console.error(`Error fetching data: ${error}`);
    res.status(500).json({ error: "Error fetching data" });
  }
}

export default handler;
