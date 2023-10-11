/* app dependencies */
import * as IconNames from "./iconNames";

/*
 Names of supported custom icons.
 */
export type IconName = (typeof IconNames)[keyof typeof IconNames];
