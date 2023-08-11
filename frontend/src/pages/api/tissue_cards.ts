import type { NextApiRequest, NextApiResponse } from "next";
import { readJson } from "src/common/utils/api/readJson";

export const allTissues = readJson(
  "src/views/CellGuide/common/fixtures/allTissues.json"
);

export const allTissueDescriptions = readJson(
  "src/views/CellGuide/common/fixtures/tissueDescriptions.json"
);

async function handler(_: NextApiRequest, res: NextApiResponse) {
  try {
    res.status(200).json(allTissues);
  } catch (error) {
    console.error(`Error fetching data: ${error}`);
    res.status(500).json({ error: "Error fetching data" });
  }
}

export default handler;
