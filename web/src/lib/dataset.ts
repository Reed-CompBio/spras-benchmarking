import { parse } from "yaml";

import scoresYaml from "../../public/data/configs/scores.yaml?raw";

const configs: Record<string, Record<string, unknown>> = {
  scores: parse(scoresYaml),
};

export const datasets = Object.entries(configs)
  .map(([type, entry]) => (entry["datasets"] as Record<string, unknown>[]).map((dataset) => ({ ...dataset, type })))
  .flat();
