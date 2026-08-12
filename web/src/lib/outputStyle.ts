interface Output {
  dataType: string;
  datasetCategory: string;
  datasetName?: string;
  algorithm: string;
  paramsHash: string;
}

function extractPrefix(name: string, prefixName: string, prefixes: string[]): { prefix: string; name: string } {
  const foundPrefix = prefixes.find((prefix) => name.startsWith(prefix));

  if (!foundPrefix) {
    throw new Error(`${name} should begin with a ${prefixName} (one of ${prefixes})`);
  }

  return {
    prefix: foundPrefix,
    name: name.substring(foundPrefix.length),
  };
}

const dataTypes = ["dmmm"];

export function extractDatasetType(name: string): { type: string; name: string } {
  const { prefix, name: newName } = extractPrefix(name, "dataset type", dataTypes);
  return { type: prefix, name: newName };
}

const dataCategories = {
  egfr: {
    name: "EGFR",
    directory: "egfr"
  },
  yeast: {
    name: "Yeast",
    directory: "yeast_osmotic_stress",
  },
};

// TODO: replace this once we have proper dataset categories
export function extractDatasetCategory(name: string): { category: string; name: string } {
  const { prefix, name: newName } = extractPrefix(name, "dataset category", Object.keys(dataCategories));
  return { category: prefix, name: newName.slice(1) };
}

export function parseOutputString(str: string): Output {
  const components = str.split("-");
  let dataType: string | undefined;
  let datasetCategory: string | undefined;
  let datasetName: string | undefined;
  let algorithm: string | undefined;
  let paramsHash: string | undefined;

  if (components.length === 5) {
    // This is a slug URL (type-...)
    [dataType, datasetCategory, datasetName, algorithm, paramsHash] = components;
  } else if (components.length === 4) {
    // This is also a slug URL w/o a name
    [dataType, datasetCategory, algorithm, paramsHash] = components;
  } else if (components.length === 3) {
    // This is fetched straight from the folder - we ignored -params previously
    [datasetName, algorithm, paramsHash] = components;
  } else {
    throw new Error(`Unexpected length of components in ${components}.`);
  }

  // We didn't get a data type in the first passthrough - lets extract the data
  // type from the name
  if (!dataType || !datasetCategory) {
    if (!datasetName) throw new Error(`datasetName ${datasetName} isn't set - this is an internal error.`);
    const { type, name: name1 } = extractDatasetType(datasetName);
    const { category, name } = extractDatasetCategory(name1);
    dataType = type;
    datasetCategory = category;
    datasetName = name;
  }

  return {
    dataType,
    datasetCategory,
    datasetName,
    algorithm,
    paramsHash,
  };
}

export function addOptional(name: string | undefined, settings: { prefix?: string; suffix?: string }): string {
  return name ? `${settings.prefix ?? ""}${name}${settings.suffix ?? ""}` : "";
}

export function styleOutput(output: Output): string {
  return `${output.dataType}-${output.datasetCategory}-${addOptional(output.datasetName, { suffix: "-" })}${output.algorithm}-${output.paramsHash}`;
}

export function asFolderName(output: Output): string {
  return `${output.dataType}${output.datasetCategory}${addOptional(output.datasetName, { prefix: "_" })}-${output.algorithm}-params-${output.paramsHash}`;
}

export function algorithmDocumentationUrl(algorithm: string): string {
  const map: Record<string, string> = {
    omicsintegrator1: "oi1",
    omicsintegrator2: "oi2",
  };

  const foundAlgorithm = algorithm in map ? map[algorithm] : algorithm;
  return `https://spras.readthedocs.io/en/latest/prms/${foundAlgorithm}.html`;
}
