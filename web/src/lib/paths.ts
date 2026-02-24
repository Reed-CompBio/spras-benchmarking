import { extractDatasetCategory, extractDatasetType } from "./outputStyle";
import { globSync } from 'glob'

export function getDataFiles(): string[] {
  // We prefer this over import.meta.glob, as import.meta.glob currently
  // leads to OOM for large raw imports, and OOM is especially plausible on CD.
  const dataFiles = globSync("../../public/data/output/**");
  return dataFiles.map((path) => path.substring("../../public/data/output/".length));
}

export function getDatasets() {
  const files = getDataFiles();
  return files
    .filter((file) => file.startsWith("logs/datasets-"))
    .map((file) => file.substring("logs/datasets-".length))
    .map((file) => file.slice(0, -".yaml".length))
    .map((file) => extractDatasetType(file))
    .map(({ type, name }) => ({ type, ...extractDatasetCategory(name) }));
}
