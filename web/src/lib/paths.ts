import { extractDatasetCategory, extractDatasetType } from "./outputStyle";
import { globSync } from 'glob'
import path from 'path';

export function getDataFiles(): string[] {
  // We prefer this over import.meta.glob, as import.meta.glob currently
  // leads to OOM for large raw imports, and OOM is especially plausible on CD.
  //(1)const dataFiles = globSync("public/data/output/**", { windowsPathsNoEscape: true });
  const dataFiles = globSync("public/data/output/**", { windowsPathsNoEscape: true });
  return dataFiles
    .map((path) => path.replace(/\\/g, "/"))
    .map((path) => path.substring("public/data/output/".length));
}


export function getDatasets() {
  const prefix=["logs", "datasets-"].join("/");
  const files = getDataFiles();
  //console.log("hello");
  //console.log(files);
  //(3)const universal= files.map((file) => file.replaceAll("\\", "/"));
  //console.log(universal);
  //console.log("bye");
  //(4)return universal
  return files
    .filter((file) => file.startsWith("logs/datasets-")) //just saves the logs\\datasets- files
    .map((file) => file.substring("logs/datasets-".length)) //[ 'dmmmyeast.yaml', 'dmmmegfr_string.yaml', 'dmmmegfr_irefindex.yaml' ]
    .map((file) => file.slice(0, -".yaml".length)) //[ 'dmmmyeast', 'dmmmegfr_string', 'dmmmegfr_irefindex' ]
    .map((file) => extractDatasetType(file)) ////[{ type: 'dmmm', name: 'yeast' },{ type: 'dmmm', name: 'egfr_string' },{ type: 'dmmm', name: 'egfr_irefindex' }]
    .map(({ type, name }) => ({ type, ...extractDatasetCategory(name) })); //[{ type: 'dmmm', category: 'yeast', name: '' },{ type: 'dmmm', category: 'egfr', name: 'string' },{ type: 'dmmm', category: 'egfr', name: 'irefindex' }]
}
