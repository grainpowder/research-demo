석사과정동안 제가 연구주제로 삼은 분야는 베이지안 통계학, 그 중에서도 변분추론(Variational Inference)이었습니다. 많은 연구들에서 이 방법을 MCMC기법과 비슷하게 정확하지만 연산 필요량을 크게 줄여 속도를 크게 향상시킬 수 있는 베이지안 추정 방법이라고 설명하고 있습니다.

2020년이 가기 전까지 완성하고자 하는 제 석사학위논문도 이 추정법에 기반한 모델을 계속해서 구현하며 준비하고 있습니다. 아직 다양한 논문을 공부하며 더 많은 기능을 구현중이기 때문에 완결되지 않은 부분이 많고, 그래서 제 연구의 모든 내용을 보여드릴 수는 없지만, 지금까지 직접 구현을 완료한 부분까지의 설명을 담기 위해 이 데모를 준비했습니다.

제가 설명드릴 모델은 *Semiparametric Mixed Effects Model with Measurement Error*라는 모델입니다. 이를 각각 Semiparametric Mixed Effects Model과 Model with Measurement Error로 나눠서 설명드리려고 합니다.

## Semiparaetric Mixed Effects Model

먼저, 혼합효과모형(Mixed Effects Model)은 주로 고정효과(Fixed Effect)와 임의효과(Random Effect)를 함께 포함하는 선형모형을 말합니다. 이 모형은 자주 환자, 기업과 같은 측정대상에 대한 반복측정데이터를 모델링하는데 사용됩니다. 이런 상황에서 임의효과는 고정효과로 설명되지 않는, 측정대상에 의한 잠재적인 효과를 추정하기 위해 포함됩니다.

이런 모형이 준모수적(Semiparametric)이라고 할 때는 이 모형에 모수적(Parametric), 비모수적(Nonparametric)요소가 혼재되어 있음을 의미합니다. 
예를 들어, 논문에 포함시킬 실제 사례연구를 위한 생체인증(Biometrics) 데이터에는 각 코호트의 소변 내 카드뮴 검출량과 특정 단백질 농도의 조화평균이 기록되어 있습니다. 두 변수 사이에 지수적으로 증가하는 양상의 패턴이 발견되는데, 이를 선형성을 가정하고 추정하는 것은 부적절하므로 이 패턴을 비모수적으로 추정하게 됩니다. 이 데이터에 포함된 다른 변수는 선형성을 가정할 수 있어서 모수적으로 추정하기 때문에, 이 모형은 준모수적입니다.

마지막으로, 임의효과에 부여되는 제약의 형태에 따라 이를 다시 두 가지로 분류할 수 있습니다.

- Random Intercept Model : 임의효과가 임의의 실수일 때
- Stochastic Frontier Model : 임의효과가 양의 실수일 때

현재 두 모형 모두 구현되어 있는 상태입니다.

## Model with Measurement Error

설명변수 x가 반응변수 y와 y = f(x) + e와 같은 관계를 갖는다고 하겠습니다(e는 랜덤 노이즈입니다). 여기에 측정오차(Measurement Error)가 포함되었다는 것은 실제 관측된 데이터가 (x, y)가 아닌 (v, y)라는 의미입니다. 예를 들어, 측정 기계가 고장이나 악화된 외부 환경 등의 이유로 정밀하지 않은 측정을 했을 때나 자료를 저장하는 과정에서 반올림을 하는 등의 이유로 정보의 손실이 일어났을 때 데이터의 참값인 x대신에 측정오차가 포함된 v를 관측하게 됩니다. 

이때 (v, y)를 통해서 (x, y)와의 관계를 파악할 수 없으므로 f의 형태에 대한 어떤 가정도 하지 않는 비모수적 추정이 적합합니다. 또한, f를 추정하는 과정에서 v의 정보를 활용해 x의 참값을 함께 추정하게 되는데, 이 과정에서 수치적분(numerical integration)이 필요하게 됩니다.

![simulation-result](./result.png)

보다 구체적인 설명을 위해 간단한 측정오차 회귀모형의 모수들을 추정하는 알고리즘을 구현해 시뮬레이션을 준비했습니다. f(x)=2x+2sin(πx)로 설정하고, y = f(x) + e라는 data generating process에 의해 y를 생성했습니다. 

이때 x에 일정량의 측정오차를 포함시켜서 생성한 v와 y의 산점도를 맨 왼쪽에 포함시켰습니다. 실제 f의 형태인 증가하는 사인파 패턴을 전혀 찾아볼 수 없습니다. 설명드리는 모델은 이처럼 측정오차로 인해 실제 패턴을 확인할 수 없게 되어버린 데이터를 가지고 숨겨진 실제 패턴을 찾아낼 때 사용할 수 있습니다.

다음 두 플랏은 추정 결과를 보여드리기 위해 포함시켰습니다. 가운데 플랏은 실제 x를 가로축에, 추정된 x를 세로축에 배치시킨 플랏입니다. y = x선이 플랏을 가로지르고 있는데, 이를 통해 추정된 x가 실제 x와 비슷하게 추정되었음을 확인할 수 있습니다. 나아가, 약 96%의 추정치들의 95% 신용구간이 실제 x를 포함하고 있는 것을 확인했고, 이를 통해서도 추정이 정확하게 이루어졌음을 확인했습니다.

마지막 플랏은 f에 대한 추정결과를 나타냅니다. 노란색 음영은 각 x에서 f(x)의 95% 신용구간, 황토색 실선은 추정된 f(x), 파란색 실선은 실제 f(x)를 의미합니다. 파란색 실선이 모든 x에서 95% 신용구간에 포함되고, 데이터가 집중된 구간에서는 증가하는 사인파 형태를 잘 추정해냈기 때문에 역시 정확한 추정이 이루어졌음을 확인했습니다. 