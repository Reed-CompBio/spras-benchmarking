import { execSync } from 'child_process';

const decoder = new TextDecoder()

export const revision = decoder.decode(execSync('git rev-parse HEAD')).trim();
export const shortRevision = revision.substring(0, 6);
