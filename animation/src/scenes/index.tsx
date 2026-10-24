import React from 'react';
import {AlleleScene} from './AlleleScene';
import {BlocksScene} from './BlocksScene';
import {CompareScene} from './CompareScene';
import {DockerScene} from './DockerScene';
import {HookScene} from './HookScene';
import {InputsScene} from './InputsScene';
import {MapScene} from './MapScene';
import {LimitsScene} from './LimitsScene';
import {OutroScene} from './OutroScene';
import {OutputsScene} from './OutputsScene';
import {ProblemScene} from './ProblemScene';
import {ProteinScene} from './ProteinScene';
import {RefCheckScene} from './RefCheckScene';
import {StopScene} from './StopScene';
import {WhyScene} from './WhyScene';
import type {SceneProps} from './types';

const SCENES: Record<string, React.FC<SceneProps>> = {
  why: WhyScene,
  map: MapScene,
  hook: HookScene,
  problem: ProblemScene,
  inputs: InputsScene,
  ref: RefCheckScene,
  allele: AlleleScene,
  blocks: BlocksScene,
  protein: ProteinScene,
  stop: StopScene,
  outputs: OutputsScene,
  docker: DockerScene,
  compare: CompareScene,
  limits: LimitsScene,
  outro: OutroScene,
};

export const SceneRenderer: React.FC<SceneProps> = (props) => {
  const Component = SCENES[props.scene.id] ?? HookScene;
  return <Component {...props} />;
};
