import { execSync } from 'child_process';

const decoder = new TextDecoder()

export const revision = decoder.decode(execSync('git rev-parse HEAD'));
export const shortRevision = revision.substring(6);
